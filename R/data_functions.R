#' Unfold EEMs into a matrix of samples and wavelength pairs
#'
#' @description Each EEM in an eemlist is unfolded into a vector. The result is a matrix with one row per sample and one column per pair of excitation and emission wavelengths. This matrix is the input for the self-organising map in \code{\link{eemsom}}. The unfolding can be reversed by \code{\link{mat2eem}}.
#'
#' @param eem_list object of class eemlist. All EEMs need to have the same excitation and emission wavelengths.
#' @param normalise normalisation of each sample: \code{"sd"} (default) divides each EEM by its standard deviation like \code{\link[staRdom]{norm_array}}, \code{"max"} by its maximum (Ejarque-Gonzalez and Butturini 2014), \code{"sum"} by the sum of its fluorescence and \code{"none"} keeps the original values.
#' @param na handling of missing values (e.g. removed scatter): \code{"remove"} (default) removes wavelength pairs that are missing in at least one sample, \code{"keep"} keeps them.
#'
#' @details Columns are named "ex|em", e.g. "250|450". The EEM matrix is unfolded column-wise, so emission wavelengths change faster than excitation wavelengths.
#'
#' Normalisation removes differences in concentration and focuses on the shape of the EEMs. The normalisation factors are calculated from the remaining wavelength pairs (after removing missing values) and are saved as an attribute.
#'
#' @return A matrix with samples in rows and wavelength pairs in columns. Attributes: \code{ex} and \code{em} (wavelengths of the original EEMs), \code{keep} (logical vector stating which of the unfolded wavelength pairs are columns of the matrix), \code{samples}, \code{normalise} and \code{norm_factors}.
#' @export
#'
#' @references Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @seealso \code{\link{mat2eem}}, \code{\link[staRdom]{eem2array}}
#'
#' @examples
#' data(eem_list)
#'
#' m <- eem2mat(eem_list)
#' dim(m)
#' attr(m, "norm_factors")
eem2mat <- function(eem_list, normalise = "sd", na = "remove"){
  eem_checkmat(eem_list)
  normalise <- match.arg(normalise, c("sd", "max", "sum", "none"))
  na <- match.arg(na, c("remove", "keep"))
  ex <- eem_list[[1]]$ex
  em <- eem_list[[1]]$em
  samples <- eem_names(eem_list)
  mat <- do.call(rbind, lapply(eem_list, function(eem) as.vector(eem$x)))
  colnames(mat) <- paste0(rep(ex, each = length(em)), "|", rep(em, times = length(ex)))
  rownames(mat) <- samples
  keep <- rep(TRUE, ncol(mat))
  if(na == "remove") keep <- !apply(is.na(mat), 2, any)
  if(!any(keep)) stop("all wavelength pairs contain missing values in at least one sample! Please consider na = \"keep\" or eem_interp().")
  mat <- mat[, keep, drop = FALSE]
  norm_factors <- eem_normfactors(mat, normalise)
  if(any(!is.finite(norm_factors) | norm_factors == 0)) stop("normalisation factors could not be calculated for the samples ", paste(samples[!is.finite(norm_factors) | norm_factors == 0], collapse = ", "), "!")
  mat <- mat / norm_factors
  attr(mat, "ex") <- ex
  attr(mat, "em") <- em
  attr(mat, "keep") <- keep
  attr(mat, "samples") <- samples
  attr(mat, "normalise") <- normalise
  attr(mat, "norm_factors") <- norm_factors
  mat
}

#' Fold a matrix of unfolded EEMs back into an eemlist
#'
#' @description Rows of a matrix created by \code{\link{eem2mat}} or of a matrix with the same columns (e.g. codebook vectors of a self-organising map) are converted to EEMs. Wavelength pairs that were removed before are set to NA.
#'
#' @param mat matrix with unfolded EEMs in rows
#' @param ex excitation wavelengths of the EEMs
#' @param em emission wavelengths of the EEMs
#' @param keep logical vector of length \code{length(ex) * length(em)} stating which wavelength pairs are columns of \code{mat}
#' @param samples sample names for the EEMs
#' @param norm_factors optional numeric vector, rows of \code{mat} are multiplied by these factors to reverse a normalisation. NULL (default) leaves the values unchanged.
#'
#' @return object of class eemlist
#' @export
#'
#' @seealso \code{\link{eem2mat}}
#'
#' @examples
#' data(eem_list)
#'
#' m <- eem2mat(eem_list)
#' eem_list2 <- mat2eem(m)                                         # normalised EEMs
#' eem_list3 <- mat2eem(m, norm_factors = attr(m, "norm_factors")) # original values
#' ggeem(eem_list2)
mat2eem <- function(mat, ex = attr(mat, "ex"), em = attr(mat, "em"), keep = attr(mat, "keep"), samples = rownames(mat), norm_factors = NULL){
  if(is.null(dim(mat))) mat <- matrix(mat, nrow = 1)
  if(is.null(ex) | is.null(em)) stop("excitation and emission wavelengths are missing!")
  if(is.null(keep)) keep <- rep(TRUE, length(ex) * length(em))
  if(length(keep) != length(ex) * length(em)) stop("length of keep does not fit the wavelengths!")
  if(sum(keep) != ncol(mat)) stop("number of columns in mat does not fit the wavelengths!")
  if(is.null(samples)) samples <- paste0("V", seq_len(nrow(mat)))
  if(!is.null(norm_factors)){
    if(length(norm_factors) != nrow(mat)) stop("length of norm_factors does not fit the number of rows in mat!")
    mat <- mat * norm_factors
  }
  lapply(seq_len(nrow(mat)), function(n){
    x <- rep(NA_real_, length(keep))
    x[keep] <- mat[n, ]
    eem <- list(sample = samples[n], x = matrix(x, nrow = length(em), ncol = length(ex)), ex = ex, em = em, location = "staRsom")
    class(eem) <- "eem"
    eem
  }) %>%
    `class<-`("eemlist")
}

#' Check whether an eemlist can be unfolded into a matrix
#'
#' @param eem_list object of class eemlist
#'
#' @return TRUE invisibly, an error is thrown otherwise
#' @noRd
eem_checkmat <- function(eem_list){
  if(!inherits(eem_list, "eemlist")) stop("eem_list is not of class eemlist!")
  if(length(eem_list) < 1) stop("eem_list is empty!")
  ex <- eem_list[[1]]$ex
  em <- eem_list[[1]]$em
  same <- vapply(eem_list, function(eem) isTRUE(all.equal(eem$ex, ex)) & isTRUE(all.equal(eem$em, em)), logical(1))
  if(!all(same)) stop("dimensions mismatch! All EEMs need the same excitation and emission wavelengths. Please consider eem_red2smallest() or eem_extend2largest().")
  invisible(TRUE)
}

#' Calculate normalisation factors for the rows of a matrix
#'
#' @param mat matrix, samples in rows
#' @param normalise normalisation method, see \code{\link{eem2mat}}
#'
#' @return numeric vector
#' @importFrom stats sd
#' @noRd
eem_normfactors <- function(mat, normalise){
  switch(normalise,
         sd = apply(mat, 1, sd, na.rm = TRUE),
         max = apply(mat, 1, max, na.rm = TRUE),
         sum = rowSums(mat, na.rm = TRUE),
         none = rep(1, nrow(mat)))
}

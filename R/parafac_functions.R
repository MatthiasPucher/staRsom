#' Calculate a self-organising map from the loadings of a PARAFAC model
#'
#' @description Instead of the EEMs, the sample loadings (A-mode) of a PARAFAC model are arranged on a self-organising map. The loadings are scaled to Fmax, the maximum fluorescence of each component in each sample. Cuss et al. (2014) found a clearer separation of DOM sources with the relative contributions of PARAFAC components than with the EEMs themselves.
#'
#' @param pfmodel object of class parafac, see \code{\link[staRdom]{eem_parafac}}
#' @param dim size of the map, either \code{"auto"} (default) or a numeric vector \code{c(xdim, ydim)}, see \code{\link{eemsom}}
#' @param normalise logical, if TRUE (default), the Fmax values of each sample are divided by their sum, so the map shows the relative composition of the samples (Cuss et al. 2014, Wheeler et al. 2017). If FALSE, the Fmax values are used directly (Zhang et al. 2015).
#' @param seed seed for the random number generator to get reproducible results
#' @param ... additional parameters passed on to \code{\link{eemsom}} (\code{topo}, \code{rlen}, \code{alpha}, \code{mode}, \code{init}, \code{cores}, \code{verbose}) and to \code{\link[kohonen]{som}}
#'
#' @details The result can be used with the functions for self-organising maps of EEMs, e.g. \code{\link{eemsom_samples}}, \code{\link{eemsom_cluster}}, \code{\link{eemsom_split}} and the plot functions. The prototype EEM of a unit (\code{\link{eemsom_codes}}, \code{\link[staRdom]{ggeem}}) is the sum of the PARAFAC components weighted by the prototype loadings. \code{\link{eemsom_plot_cplane}} shows component planes of the PARAFAC components with the argument \code{comps}. \code{\link{eemsom_predict}} maps new samples after calculating their loadings with the fixed emission and excitation loadings of the PARAFAC model.
#'
#' Relative contributions (\code{normalise = TRUE}) are only meaningful for non-negative loadings, i.e. for models calculated with \code{const = c("nonneg", "nonneg", "nonneg")} in \code{\link[staRdom]{eem_parafac}}.
#'
#' \code{\link{eemsom_regions}} and \code{\link{eemsom_stability}} are not available for these maps.
#'
#' @return object of class eempf_som and eemsom, see \code{\link{eemsom}}. Additionally, it contains the element pfmodel, the PARAFAC model with emission (B) and excitation (C) loadings scaled to Fmax (see \code{\link[staRdom]{eempf_rescaleBC}}), which is used to reconstruct prototype EEMs and to map new samples. The element normalise is the logical value of the argument and norm_factors are the sums of the Fmax values (or 1).
#' @export
#'
#' @importFrom staRdom eempf_rescaleBC norm2A A_missing
#'
#' @references Cuss, C. W., Shi, Y. X., McConnell, S. M., Gueguen, C. (2014). Changes in the fluorescence composition of multiple DOM sources over pH gradients assessed by combining parallel factor analysis and self-organizing maps. Journal of Geophysical Research: Biogeosciences 119(9), 1850-1860. \doi{10.1002/2013JG002598}
#'
#' Wheeler, K. I., Levia, D. F., Hudson, J. E. (2017). Tracking senescence-induced patterns in leaf litter leachate using parallel factor analysis (PARAFAC) modeling and self-organizing maps. Journal of Geophysical Research: Biogeosciences 122(9), 2233-2250. \doi{10.1002/2016JG003677}
#'
#' Zhang, Y., Liang, X., Wang, Z., Xu, L. (2015). A novel approach combining self-organizing map and parallel factor analysis for monitoring water quality of watersheds under non-point source pollution. Scientific Reports 5, 16079. \doi{10.1038/srep16079}
#'
#' @seealso \code{\link{eemsom}}, \code{\link[staRdom]{eem_parafac}}
#' @family self-organising maps
#'
#' @examples
#' data(pf_models)
#'
#' pfsom <- eempf_som(pf4[[1]], rlen = 100)
#' pfsom
#' eemsom_plot_cplane(pfsom, comps = 1:6)
#' pfsom <- eemsom_cluster(pfsom, k = 4)
#' eemsom_plot_map(pfsom, type = "cluster")
#' ggeem(pfsom, units = 1:4)
eempf_som <- function(pfmodel, dim = "auto", normalise = TRUE, seed = 42, ...){
  if(!inherits(pfmodel, "parafac")) stop("pfmodel must be an object of class parafac!")
  # B and C scaled to Fmax, the normalisation of the EEMs is kept in the model for A_missing()
  pf <- eempf_rescaleBC(pfmodel)
  loadings <- norm2A(pf)$A
  if(is.null(rownames(loadings))) stop("The loadings of pfmodel have no sample names!")
  scaled <- eempf_normloadings(loadings, normalise)
  trained <- eemsom_train(scaled$mat, dim = dim, seed = seed, varname = "components", ...)
  res <- list(som = trained$som, ex = as.numeric(rownames(pf$C)), em = as.numeric(rownames(pf$B)), keep = rep(TRUE, nrow(pf$C) * nrow(pf$B)),
              samples = rownames(loadings), normalise = normalise, norm_factors = scaled$norm_factors, na = "remove", clusters = NULL,
              settings = trained$settings, pfmodel = pf)
  class(res) <- c("eempf_som", "eemsom")
  res
}

#' Scale Fmax loadings for a self-organising map
#'
#' @param loadings matrix of Fmax loadings, samples in rows
#' @param normalise logical, see eempf_som
#'
#' @return list with mat (relative or absolute loadings) and norm_factors
#' @noRd
eempf_normloadings <- function(loadings, normalise){
  norm_factors <- if(normalise) rowSums(loadings) else rep(1, nrow(loadings))
  if(any(norm_factors <= 0)) stop("The loadings of the following samples cannot be normalised (sum <= 0): ", paste(rownames(loadings)[norm_factors <= 0], collapse = ", "))
  list(mat = loadings / norm_factors, norm_factors = unname(norm_factors))
}

#' Map new samples to a self-organising map of PARAFAC loadings
#'
#' @param som object of class eempf_som
#' @param eem_list object of class eemlist with the wavelengths of the PARAFAC model
#' @param cores number of cores for A_missing
#' @param ... passed on to A_missing
#'
#' @return data frame, see eemsom_predict
#' @noRd
eempf_som_predict <- function(som, eem_list, cores = 1, ...){
  if(is.null(som$pfmodel)) stop("som does not contain its PARAFAC model, please calculate it again with eempf_som()!")
  eem_checkmat(eem_list)
  same_wl <- isTRUE(all.equal(eem_list[[1]]$ex, som$ex)) && isTRUE(all.equal(eem_list[[1]]$em, som$em))
  if(!same_wl) stop("The wavelengths of eem_list differ from the PARAFAC model! Please consider eem_range(), eem_red2smallest() or eem_extend2largest().")
  # the loadings with fixed B and C modes do not depend on the start, so one start is enough
  args <- utils::modifyList(list(eem_list = eem_list, pfmodel = som$pfmodel, cores = cores, nstart = 1), list(...))
  # A_missing() warns about surplus wavelengths even if there are none and prints a progress bar, both are suppressed
  invisible(utils::capture.output(
    pf <- withCallingHandlers(do.call(A_missing, args), warning = function(w){
      if(grepl("not present in the PARAFAC model", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  ))
  loadings <- pf$A
  # loadings of the original (not normalised) EEMs, like norm2A()
  if(!is.null(attr(pf, "norm_factors"))) loadings <- loadings * attr(pf, "norm_factors")
  mat <- eempf_normloadings(loadings, som$normalise)$mat
  d <- eemsom_distmat(mat, som$som$codes[[1]])
  unit <- apply(d, 1, which.min)
  eemsom_unittable(som, rownames(mat), unit, d[cbind(seq_along(unit), unit)])
}

#' @rdname eempf_som
#' @param x object of class eempf_som
#' @export
print.eempf_som <- function(x, ...){
  cat("Self-organising map of", length(x$samples), "samples from a PARAFAC model with", ncol(x$pfmodel$B), "components\n")
  cat(" map:", x$som$grid$xdim, "x", x$som$grid$ydim, x$som$grid$topo, "units\n")
  cat(" loadings:", if(x$normalise) "relative (Fmax / sum of Fmax)" else "Fmax", "\n")
  cat(" clusters:", if(is.null(x$clusters)) "none" else length(unique(x$clusters)), "\n")
  invisible(x)
}

#' Stop if a function is not available for maps of PARAFAC loadings
#'
#' @param som object of class eemsom
#' @param fun name of the function
#'
#' @return TRUE invisibly
#' @noRd
eempf_som_unsupported <- function(som, fun){
  if(inherits(som, "eempf_som")) stop(fun, "() is not available for self-organising maps of PARAFAC loadings!")
  invisible(TRUE)
}

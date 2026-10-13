#' Calculate a self-organising map from EEM data
#'
#' @description The EEMs are unfolded into vectors (see \code{\link{eem2mat}}), normalised and arranged on a self-organising map (SOM, Kohonen 2013). Each unit of the map gets a prototype EEM (codebook vector). Similar samples are mapped to the same or neighbouring units. The calculation is done with \code{\link[kohonen]{som}}, please see details there.
#'
#' @param eem_list object of class eemlist. All EEMs need to have the same excitation and emission wavelengths and unique sample names.
#' @param dim size of the map, either \code{"auto"} (default) or a numeric vector \code{c(xdim, ydim)} with the number of columns and rows of units
#' @param topo topology of the map, \code{"hexagonal"} (default) or \code{"rectangular"}
#' @param normalise normalisation of each sample, see \code{\link{eem2mat}}
#' @param na handling of missing values, see \code{\link{eem2mat}}
#' @param rlen number of times the complete data set is presented to the map during training
#' @param alpha learning rate, vector of start and end value
#' @param mode training algorithm: \code{"online"} (default), \code{"batch"} or \code{"pbatch"} (parallel batch)
#' @param init initialisation of the codebook vectors: \code{"pca"} (default) places them on the plane of the first two principal components, \code{"random"} uses randomly selected samples
#' @param seed seed for the random number generator to get reproducible results. NULL does not set a seed. The random number generator state of the session is restored afterwards.
#' @param cores number of cores used with \code{mode = "pbatch"}
#' @param verbose print infos
#' @param ... additional parameters that are passed on to \code{\link[kohonen]{som}}, e.g. \code{radius}. \code{keep.data} is always TRUE, because the data are needed by \code{\link{eemsom_samples}} and other functions.
#'
#' @details With \code{dim = "auto"}, the number of units is \eqn{5 \sqrt{n}}{5 sqrt(n)} (Vesanto and Alhoniemi 2000) but not more than the number of samples n. The ratio of the side lengths follows the two largest eigenvalues of the data (Bieroza et al. 2009, Ejarque-Gonzalez and Butturini 2014): like in the SOM Toolbox (Vesanto et al. 2000), it is the square root of the ratio of the eigenvalues, i.e. the ratio of the standard deviations along the first two principal components. Ratios above the square root of the number of units are replaced by 1 to avoid very elongated maps.
#'
#' The default learning rate (0.05 to 0.01) and the default neighbourhood radius (decreasing from two thirds of all unit distances to zero) of \code{\link[kohonen]{som}} correspond to the settings of Ejarque-Gonzalez and Butturini (2014).
#'
#' The squared Euclidean distance (\code{dist.fcts = "sumofsquares"}) is always used. Otherwise, kohonen would switch to the Tanimoto distance for data between 0 and 1 with row sums of 1, which occurs with \code{normalise = "sum"}. With \code{na = "keep"}, distances are calculated from the available wavelength pairs of a sample and scaled up to the number of all wavelength pairs.
#'
#' The initialisation with \code{init = "pca"} is deterministic. The order of samples in the online training is random, so \code{seed} is needed for reproducible results.
#'
#' The training time grows with the numbers of samples, units and wavelength pairs. As a rough guide for one core, online training with \code{rlen = 500} takes about a minute for 200 EEMs with 3000 wavelength pairs and several minutes for 1000 EEMs. \code{mode = "pbatch"} with \code{cores > 1}, a smaller \code{rlen} or a lower resolution of the EEMs (see \code{\link[staRdom]{eem_range}} and \code{\link[staRdom]{eem_exclude}}; Ejarque-Gonzalez and Butturini (2014) used 366 wavelength pairs) reduce the time. This matters for \code{\link{eemsom_stability}}, which repeats the training many times.
#'
#' @return object of class eemsom, a list containing
#' \itemize{
#'   \item som: object of class kohonen, see \code{\link[kohonen]{som}}
#'   \item ex, em: excitation and emission wavelengths
#'   \item keep: logical vector of the unfolded wavelength pairs used in the SOM
#'   \item samples: sample names
#'   \item normalise, norm_factors: normalisation method and factors
#'   \item na: handling of missing values
#'   \item clusters: cluster of each unit, NULL until the units are clustered
#'   \item settings: arguments used for \code{\link[kohonen]{som}}, needed by \code{\link{eemsom_stability}}
#' }
#' @export
#'
#' @import kohonen
#'
#' @references Kohonen, T. (2013). Essentials of the self-organizing map. Neural Networks 37, 52-65. \doi{10.1016/j.neunet.2012.09.018}
#'
#' Vesanto, J., Alhoniemi, E. (2000). Clustering of the self-organizing map. IEEE Transactions on Neural Networks 11(3), 586-600. \doi{10.1109/72.846731}
#'
#' Vesanto, J., Himberg, J., Alhoniemi, E., Parhankangas, J. (2000). SOM Toolbox for Matlab 5. Report A57, Helsinki University of Technology.
#'
#' Bieroza, M., Baker, A., Bridgeman, J. (2009). Exploratory analysis of excitation-emission matrix fluorescence spectra with self-organizing maps as a basis for determination of organic matter removal efficiency at water treatment works. Journal of Geophysical Research: Biogeosciences 114(G4), 2009JG000940. \doi{10.1029/2009JG000940}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' Wehrens, R., Kruisselbrink, J. (2018). Flexible Self-Organizing Maps in kohonen 3.0. Journal of Statistical Software 87(7), 1-18. \doi{10.18637/jss.v087.i07}
#'
#' @seealso \code{\link[kohonen]{som}}, \code{\link{eemsom_samples}}, \code{\link{eemsom_codes}}, \code{\link{eemsom_quality}}
#' @family self-organising maps
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100, seed = 42)
#' som
#' eemsom_samples(som)
eemsom <- function(eem_list, dim = "auto", topo = "hexagonal", normalise = "sd", na = "remove", rlen = 500, alpha = c(0.05, 0.01), mode = "online", init = "pca", seed = 42, cores = 1, verbose = FALSE, ...){
  na <- match.arg(na, c("remove", "keep"))
  mat <- eem2mat(eem_list, normalise = normalise, na = na)
  extra <- list(...)
  if(na == "keep" && is.null(extra$maxNA.fraction)) extra$maxNA.fraction <- 1 - 1 / ncol(mat)
  trained <- do.call(eemsom_train, c(list(mat = mat, dim = dim, topo = topo, rlen = rlen, alpha = alpha, mode = mode, init = init, seed = seed, cores = cores, verbose = verbose, varname = "wavelength pairs"), extra))
  eemsom <- list(som = trained$som, ex = attr(mat, "ex"), em = attr(mat, "em"), keep = attr(mat, "keep"), samples = attr(mat, "samples"),
                 normalise = attr(mat, "normalise"), norm_factors = attr(mat, "norm_factors"), na = na, clusters = NULL,
                 settings = trained$settings)
  class(eemsom) <- "eemsom"
  eemsom
}

#' Train a self-organising map on a matrix
#'
#' @param mat matrix, samples in rows
#' @param dim,topo,rlen,alpha,mode,init,seed,cores,verbose see eemsom
#' @param varname name of the variables for verbose output
#' @param ... passed on to kohonen::som
#'
#' @return list with som (kohonen object) and settings
#' @noRd
eemsom_train <- function(mat, dim = "auto", topo = "hexagonal", rlen = 500, alpha = c(0.05, 0.01), mode = "online", init = "pca", seed = 42, cores = 1, verbose = FALSE, varname = "variables", ...){
  topo <- match.arg(topo, c("hexagonal", "rectangular"))
  mode <- match.arg(mode, c("online", "batch", "pbatch"))
  init <- match.arg(init, c("pca", "random"))
  if(nrow(mat) < 3) stop("at least 3 samples are needed for a self-organising map!")
  # the principal components are needed for the map size and for the initialisation
  pca <- if(identical(dim, "auto") || init == "pca") eemsom_pca(mat) else NULL
  if(identical(dim, "auto")){
    dim <- eemsom_dim(pca, nrow(mat), topo = topo)
  } else if(!is.numeric(dim) || length(dim) != 2 || any(dim < 1)){
    stop("dim needs to be \"auto\" or a numeric vector c(xdim, ydim)!")
  }
  if(prod(dim) > nrow(mat)) warning("The map has more units (", prod(dim), ") than samples (", nrow(mat), ")!")
  if(verbose) cat(paste0("A map of ", dim[1], " x ", dim[2], " units is calculated from ", nrow(mat), " samples and ", ncol(mat), " ", varname, "."), fill = TRUE)
  grid <- somgrid(xdim = dim[1], ydim = dim[2], topo = topo)
  # kohonen uses the Tanimoto distance for data in [0, 1] with row sums of 1 (e.g. normalise = "sum" or relative PARAFAC loadings), so the distance is set explicitly
  args <- utils::modifyList(list(X = mat, grid = grid, rlen = rlen, alpha = alpha, mode = mode, cores = cores, dist.fcts = "sumofsquares"), list(...))
  # the data are needed by eemsom_samples(), eemsom_quality() and the plot functions
  args$keep.data <- TRUE
  settings <- c(args[setdiff(names(args), c("X", "init", "keep.data"))], list(init = init))
  if(init == "pca") args$init <- eemsom_pcainit(pca, grid, colnames(mat))
  som <- with_seed(seed, do.call(kohonen::som, args))
  list(som = som, settings = settings)
}

#' @rdname eemsom
#' @param x object of class eemsom
#' @export
print.eemsom <- function(x, ...){
  cat("Self-organising map of", length(x$samples), "EEMs\n")
  cat(" map:", x$som$grid$xdim, "x", x$som$grid$ydim, x$som$grid$topo, "units\n")
  cat(" wavelength pairs:", sum(x$keep), "of", length(x$keep), "\n")
  cat(" normalisation:", x$normalise, "\n")
  cat(" clusters:", if(is.null(x$clusters)) "none" else length(unique(x$clusters)), "\n")
  invisible(x)
}

#' Assignment of samples to units of a self-organising map
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#'
#' @details With missing values (\code{na = "keep"}), the distances are calculated from the available wavelength pairs of a sample and scaled up to the number of all wavelength pairs, as in kohonen.
#'
#' @details The columns hits and umatrix help to find unusual samples (outliers): such samples often occupy a unit alone (hits = 1) that is separated from its neighbours by a high U-matrix value. Their distance to the prototype can be small though, because the prototype adapts to them. The leave-one-out distance calculated by \code{\link{eemsom_stability}} with \code{type = "loo"} is a more reliable indicator.
#'
#' @return data frame with columns sample, unit (number of the best matching unit), x and y (position of the unit on the map), distance (Euclidean distance between sample and prototype), hits (number of samples mapped to the unit), umatrix (mean distance between the prototype of the unit and the prototypes of its neighbours, see \code{\link{eemsom_plot_map}}) and cluster (NA if no clustering was done)
#' @export
#'
#' @family map results
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' eemsom_samples(som)
eemsom_samples <- function(som){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  eemsom_unittable(som, som$samples, som$som$unit.classif, som$som$distances)
}

#' Table of samples and their units
#'
#' @param som object of class eemsom
#' @param samples sample names
#' @param unit best matching unit of each sample
#' @param sqdist squared distance between each sample and its prototype
#'
#' @return data frame, see eemsom_samples
#' @noRd
eemsom_unittable <- function(som, samples, unit, sqdist){
  pts <- som$som$grid$pts
  hits <- tabulate(som$som$unit.classif, nbins = nrow(pts))
  cluster <- if(is.null(som$clusters)) NA else som$clusters[unit]
  data.frame(sample = samples, unit = unit, x = pts[unit, 1], y = pts[unit, 2], distance = sqrt(sqdist),
             hits = hits[unit], umatrix = eemsom_umatrix(som)[unit], cluster = cluster, stringsAsFactors = FALSE)
}

#' Prototype EEMs of a self-organising map
#'
#' @description The codebook vectors of the units are converted to EEMs.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param denormalise logical, if TRUE, each prototype is multiplied by the mean normalisation factor of the samples mapped to the unit (the mean of all samples for empty units). The result has the scale of the original data.
#'
#' @details For maps of PARAFAC loadings (\code{\link{eempf_som}}), the prototype EEM of a unit is the sum of the PARAFAC components weighted by the prototype loadings.
#'
#' @return object of class eemlist, samples are named by the unit (e.g. "V1"), see \code{\link{mat2eem}}
#' @export
#'
#' @family map results
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' codes <- eemsom_codes(som)
#' codes
#' ggeem(codes)
eemsom_codes <- function(som, denormalise = FALSE){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  codes <- som$som$codes[[1]]
  norm_factors <- NULL
  if(denormalise){
    unit <- som$som$unit.classif
    norm_factors <- sapply(seq_len(nrow(codes)), function(u){
      if(any(unit == u)) mean(som$norm_factors[unit == u]) else mean(som$norm_factors)
    })
  }
  if(inherits(som, "eempf_som")){
    # prototype EEM = sum of PARAFAC components weighted by the prototype loadings
    B <- som$pfmodel$B
    C <- som$pfmodel$C
    codes <- t(apply(codes, 1, function(a) as.vector(B %*% (a * t(C))))) %>%
      `rownames<-`(rownames(codes))
  }
  mat2eem(codes, ex = som$ex, em = som$em, keep = som$keep, samples = rownames(codes), norm_factors = norm_factors)
}

#' Quality measures of a self-organising map
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#'
#' @details
#' \itemize{
#'   \item quantisation_error: mean Euclidean distance between samples and their prototypes, see \code{\link{eemsom_samples}} for the handling of missing values
#'   \item topographic_error: share of samples whose best and second best matching units are not neighbours (six neighbours on hexagonal maps, eight including the diagonals on rectangular maps, see \code{\link[kohonen]{unit.distances}})
#'   \item explained_variance: share of the total sum of squares of the (normalised) data explained by the prototypes
#' }
#'
#' @return data frame with one row
#' @export
#'
#' @references Kohonen, T. (2013). Essentials of the self-organizing map. Neural Networks 37, 52-65. \doi{10.1016/j.neunet.2012.09.018}
#'
#' @family map results
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' eemsom_quality(som)
eemsom_quality <- function(som){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  mat <- som$som$data[[1]]
  codes <- som$som$codes[[1]]
  d <- eemsom_distmat(mat, codes)
  dmin <- apply(d, 1, min)
  bmu2 <- t(apply(d, 1, function(r) order(r)[1:2]))
  te <- mean(eemsom_unitdist(som$som$grid)[bmu2] > 1)
  # total sum of squares, scaled for missing values like the distances
  available <- !is.na(mat)
  ss_rows <- rowSums(sweep(mat, 2, colMeans(mat, na.rm = TRUE))^2, na.rm = TRUE) * ncol(mat) / rowSums(available)
  data.frame(quantisation_error = mean(sqrt(dmin)), topographic_error = te, explained_variance = 1 - sum(dmin) / sum(ss_rows))
}

#' Squared Euclidean distances between samples and codebook vectors
#'
#' @param mat matrix, samples in rows, may contain NA
#' @param codes matrix, codebook vectors in rows
#'
#' @return matrix samples x units. Like in kohonen, the distances of samples with missing values are calculated from the available variables and scaled up to the number of all variables.
#' @noRd
eemsom_distmat <- function(mat, codes){
  available <- !is.na(mat)
  x <- mat
  x[!available] <- 0
  # sum((x - c)^2) over the available variables, expanded into matrix products
  d <- rowSums(x^2) - 2 * x %*% t(codes) + available %*% t(codes^2)
  d[d < 0] <- 0
  unname(d * (ncol(mat) / rowSums(available)))
}

#' Distances between the units of a map
#'
#' @param grid object of class somgrid
#'
#' @return matrix, distances are rounded because neighbours on hexagonal maps are only approximately 1 apart
#' @noRd
eemsom_unitdist <- function(grid){
  round(unit.distances(grid), 6)
}

#' Map size according to Vesanto and Alhoniemi (2000)
#'
#' @param pca principal components of the data, see eemsom_pca
#' @param n number of samples
#' @param topo topology
#'
#' @return numeric vector c(xdim, ydim)
#' @noRd
eemsom_dim <- function(pca, n, topo = "hexagonal"){
  munits <- min(round(5 * sqrt(n)), n)
  d <- pca$d
  # ratio of the side lengths = ratio of the standard deviations along the first two principal components,
  # extreme ratios are replaced by 1 like in the SOM Toolbox (som_topol_struct)
  ratio <- if(length(d) > 1 && d[2] > 0 && d[1] <= sqrt(munits) * d[2]) d[1] / d[2] else 1
  xdim <- max(2, round(sqrt(munits / ratio * if(topo == "hexagonal") sqrt(0.75) else 1)))
  ydim <- max(2, round(munits / xdim))
  # rounding must not lead to more units than samples
  if(xdim * ydim > n) ydim <- max(2, floor(n / xdim))
  c(xdim, ydim)
}

#' Principal components of a matrix with missing values replaced by column means
#'
#' @param mat matrix, samples in rows
#'
#' @return list with centre, d (singular values) and v (loadings) of the first two components and n (number of samples)
#' @noRd
eemsom_pca <- function(mat){
  centre <- colMeans(mat, na.rm = TRUE)
  mc <- sweep(mat, 2, centre)
  mc[is.na(mc)] <- 0
  sv <- svd(mc, nu = 0, nv = min(2, dim(mc)))
  list(centre = centre, d = sv$d[seq_len(ncol(sv$v))], v = sv$v, n = nrow(mat))
}

#' Linear initialisation of codebook vectors on the plane of the first two principal components
#'
#' @param pca principal components of the data, see eemsom_pca
#' @param grid object of class somgrid
#' @param varnames column names of the result
#'
#' @return matrix of initial codebook vectors
#' @noRd
eemsom_pcainit <- function(pca, grid, varnames = NULL){
  pts <- grid$pts
  # unit coordinates scaled to [-1, 1]
  coords <- apply(pts, 2, function(p) if(diff(range(p)) > 0) 2 * (p - min(p)) / diff(range(p)) - 1 else p * 0)
  # the longer side of the map follows the first principal component
  if(diff(range(pts[, 2])) > diff(range(pts[, 1]))) coords <- coords[, 2:1]
  k <- ncol(pca$v)
  init <- matrix(pca$centre, nrow = nrow(pts), ncol = length(pca$centre), byrow = TRUE) +
    coords[, seq_len(k), drop = FALSE] %*% (t(pca$v) * pca$d / sqrt(pca$n - 1))
  colnames(init) <- varnames
  init
}

#' Evaluate code with a seed and restore the random number generator state afterwards
#'
#' @param seed integer or NULL
#' @param code code to evaluate
#'
#' @return result of code
#' @noRd
with_seed <- function(seed, code){
  if(is.null(seed)) return(code)
  # .Random.seed lives in the global environment and is restored like withr::with_seed does
  genv <- globalenv()
  old_seed <- if(exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL
  on.exit({
    if(is.null(old_seed)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old_seed, envir = genv)
  })
  set.seed(seed)
  code
}

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
#' @details The result can be used with the functions for self-organising maps of EEMs, e.g. \code{\link{eemsom_samples}}, \code{\link{eemsom_cluster}}, \code{\link{eemsom_split}} and the plot functions. The prototype EEM of a unit (\code{\link{eemsom_codes}}, \code{\link[staRdom]{ggeem}}) is the sum of the PARAFAC components weighted by the prototype loadings. \code{\link{eemsom_plot_cplane}} shows component planes of the PARAFAC components with the argument \code{comps}.
#'
#' Relative contributions (\code{normalise = TRUE}) are only meaningful for non-negative loadings, i.e. for models calculated with \code{const = c("nonneg", "nonneg", "nonneg")} in \code{\link[staRdom]{eem_parafac}}.
#'
#' \code{\link{eemsom_regions}}, \code{\link{eemsom_predict}} and \code{\link{eemsom_stability}} are not available for these maps.
#'
#' @return object of class eempf_som and eemsom, see \code{\link{eemsom}}. Additionally, it contains the element components with the Fmax-scaled emission (B) and excitation (C) loadings of the PARAFAC model. The element normalise is the logical value of the argument and norm_factors are the sums of the Fmax values (or 1).
#' @export
#'
#' @importFrom staRdom eempf_rescaleBC norm2A
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
  pf <- pfmodel %>%
    eempf_rescaleBC() %>%
    norm2A()
  loadings <- pf$A
  if(is.null(rownames(loadings))) stop("The loadings of pfmodel have no sample names!")
  norm_factors <- if(normalise) rowSums(loadings) else rep(1, nrow(loadings))
  if(any(norm_factors <= 0)) stop("The loadings of the following samples cannot be normalised (sum <= 0): ", paste(rownames(loadings)[norm_factors <= 0], collapse = ", "))
  mat <- loadings / norm_factors
  trained <- eemsom_train(mat, dim = dim, seed = seed, varname = "components", ...)
  res <- list(som = trained$som, ex = as.numeric(rownames(pf$C)), em = as.numeric(rownames(pf$B)), keep = rep(TRUE, nrow(pf$C) * nrow(pf$B)),
              samples = rownames(loadings), normalise = normalise, norm_factors = unname(norm_factors), na = "remove", clusters = NULL,
              settings = trained$settings, components = list(B = pf$B, C = pf$C))
  class(res) <- c("eempf_som", "eemsom")
  res
}

#' @rdname eempf_som
#' @param x object of class eempf_som
#' @export
print.eempf_som <- function(x, ...){
  cat("Self-organising map of", length(x$samples), "samples from a PARAFAC model with", ncol(x$components$B), "components\n")
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

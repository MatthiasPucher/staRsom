#' Find fluorescence regions from correlated component planes
#'
#' @description Wavelength pairs with similar behaviour across the units of a self-organising map are grouped into regions. Each region is assumed to represent a fluorescence component (Ejarque-Gonzalez and Butturini 2014). The component planes are arranged on a second self-organising map ("R-mode") and its units are clustered hierarchically.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param k number of regions (a single number) or candidates (a vector, default 2 to 9) from which the number with the highest mean silhouette width is chosen. NULL is the same as \code{2:9}.
#' @param method \code{"planes"} (default) or \code{"cormat"}, see details
#' @param dim size of the R-mode map, either \code{"auto"} (default) or a numeric vector \code{c(xdim, ydim)}, see \code{\link{eemsom}}
#' @param cluster_method clustering method of the R-mode units, see \code{\link{eemsom_cluster}}
#' @param seed seed for the random number generator to get reproducible results
#' @param ... additional parameters that are passed on to \code{\link[kohonen]{som}}. The default number of iterations is \code{rlen = 100}.
#'
#' @details Following Barreto-Sanz and Perez-Uribe (2007), Ejarque-Gonzalez and Butturini (2014) calculated a self-organising map from the matrix of Pearson's correlations between all component planes (\code{method = "cormat"}). This matrix has as many rows and columns as wavelength pairs and the calculation gets slow with a high resolution of the EEMs (minutes with a few thousand wavelength pairs).
#'
#' With \code{method = "planes"}, each component plane is centred and scaled to length 1. The squared Euclidean distance of two such vectors equals \eqn{2(1 - r)}{2(1 - r)}, where r is their correlation. The R-mode map is calculated from these vectors directly, which is much faster and groups wavelength pairs by their correlation as well.
#'
#' Both methods group wavelength pairs by the correlation of their component planes, but the R-mode maps are trained on different vectors. The resulting regions are similar but not identical and the region numbers can differ.
#'
#' Wavelength pairs with constant component planes have no correlation and get no region (NA), as well as wavelength pairs that were not used in the map.
#'
#' @return object of class eemsom_regions, a list containing
#' \itemize{
#'   \item regions: data frame with columns ex, em and region
#'   \item k: number of regions
#'   \item kselect: data frame with criteria for the candidates (see \code{\link{eemsom_kselect}}), NULL if k was a single number
#'   \item rsom: R-mode self-organising map of class kohonen
#'   \item unit_regions: region of each unit of the R-mode map
#'   \item method, cluster_method, seed: settings used
#' }
#' @export
#'
#' @references Barreto-Sanz, M. A., Perez-Uribe, A. (2007). Improving the Correlation Hunting in a Large Quantity of SOM Component Planes. In: Artificial Neural Networks - ICANN 2007, Lecture Notes in Computer Science 4669, 379-388. \doi{10.1007/978-3-540-74695-9_39}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @seealso \code{\link{eemsom_fri}}
#' @family fluorescence regions
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' regions <- eemsom_regions(som, dim = c(8, 8))
#' regions
#' regions$kselect
#' ggeem(regions)
#' # fixed number of regions
#' eemsom_regions(som, k = 3, dim = c(8, 8))
eemsom_regions <- function(som, k = 2:9, method = "planes", dim = "auto", cluster_method = "complete", seed = 42, ...){
  if(!inherits(som, "eemsom")) stop("som needs to be of class eemsom!")
  eempf_som_unsupported(som, "eemsom_regions")
  method <- match.arg(method, c("planes", "cormat"))
  if(is.null(k)) k <- 2:9
  if(!is.numeric(k) || length(k) < 1 || any(k < 1)) stop("k needs to be a number of regions or a vector of candidates!")
  codes <- som$som$codes[[1]]
  if(nrow(codes) < 3) stop("The self-organising map needs at least 3 units!")
  planes <- t(codes)
  planes <- planes - rowMeans(planes)
  len <- sqrt(rowSums(planes^2))
  # constant component planes (up to rounding errors) have no correlation
  valid <- len > 1e-8 * max(len)
  planes <- planes[valid, , drop = FALSE] / len[valid]
  if(method == "cormat"){
    if(nrow(planes) > 1500) warning("method = \"cormat\" with ", nrow(planes), " wavelength pairs is slow. Please consider method = \"planes\" or a lower resolution of the EEMs.")
    planes <- tcrossprod(planes)
  }
  pca <- eemsom_pca(planes)
  if(identical(dim, "auto")){
    dim <- eemsom_dim(pca, nrow(planes), topo = "hexagonal")
  } else if(!is.numeric(dim) || length(dim) != 2 || any(dim < 1)){
    stop("dim needs to be \"auto\" or a numeric vector c(xdim, ydim)!")
  }
  grid <- somgrid(xdim = dim[1], ydim = dim[2], topo = "hexagonal")
  args <- utils::modifyList(list(X = planes, grid = grid, rlen = 100, init = eemsom_pcainit(pca, grid, colnames(planes)), dist.fcts = "sumofsquares"), list(...))
  rsom <- with_seed(seed, do.call(kohonen::som, args))
  rcodes <- rsom$codes[[1]]
  kselect <- NULL
  if(length(k) > 1){
    kselect <- eemsom_kcriteria(rcodes, k = k, method = cluster_method, seed = seed)
    k <- kselect$k[which.max(kselect$silhouette)]
  } else if(k > nrow(rcodes)){
    stop("k needs to be a number between 1 and the number of units of the R-mode map (", nrow(rcodes), ")!")
  }
  unit_regions <- eemsom_clusterunits(rcodes, k = k, method = cluster_method, seed = seed)
  region <- rep(NA_integer_, length(som$keep))
  region[which(som$keep)[valid]] <- unit_regions[rsom$unit.classif]
  regions <- data.frame(ex = rep(som$ex, each = length(som$em)), em = rep(som$em, times = length(som$ex)), region = region)
  res <- list(regions = regions, k = k, kselect = kselect, rsom = rsom, unit_regions = unit_regions, method = method, cluster_method = cluster_method, seed = seed)
  class(res) <- "eemsom_regions"
  res
}

#' @rdname eemsom_regions
#' @param x object of class eemsom_regions
#' @export
print.eemsom_regions <- function(x, ...){
  cat("Fluorescence regions from a self-organising map\n")
  cat(" regions:", x$k, if(is.null(x$kselect)) "(given)" else "(highest mean silhouette)", "\n")
  cat(" wavelength pairs per region:\n")
  print(table(region = x$regions$region, useNA = "ifany"))
  invisible(x)
}

#' Plot fluorescence regions
#'
#' @description The regions found by \code{\link{eemsom_regions}} are shown in the excitation-emission space.
#'
#' @param data object of class eemsom_regions
#' @param fill_max not used, for consistency with \code{\link[staRdom]{ggeem}}
#' @param ... further arguments, not used
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @family fluorescence regions
#' @family plots
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' regions <- eemsom_regions(som, dim = c(8, 8))
#' ggeem(regions)
ggeem.eemsom_regions <- function(data, fill_max = FALSE, ...){
  data$regions %>%
    mutate(region = factor(region)) %>%
    ggplot(aes(x = ex, y = em, fill = region)) +
    geom_tile() +
    scale_fill_discrete(na.value = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Excitation (nm)", y = "Emission (nm)")
}

#' Fluorescence regional integration
#'
#' @description The fluorescence of each sample is integrated within each region. Following Chen et al. (2003), the volumes are normalised by the area of the region, so regions of different size can be compared. Ejarque-Gonzalez and Butturini (2014) quantified their fluorescence components from the self-organising map in this way.
#'
#' @param eem_list object of class eemlist with the same wavelengths as the self-organising map used for the regions. Usually, these are the original (not normalised) EEMs.
#' @param regions object of class eemsom_regions, see \code{\link{eemsom_regions}}
#'
#' @details The volume of region i is \eqn{\Phi_i = \sum I(ex, em) \Delta ex \Delta em}{Phi_i = sum(I(ex, em) * dex * dem)} over all wavelength pairs of the region. The area-normalised volume is \eqn{\Phi_i \cdot A / A_i}{Phi_i * A / A_i}, with A being the total area of all regions and \eqn{A_i}{A_i} the area of region i. Missing values in a sample are ignored and the area of the region is reduced accordingly.
#'
#' @return data frame with column sample and one column per region (region_1, region_2, ...) containing the area-normalised volumes
#' @export
#'
#' @references Chen, W., Westerhoff, P., Leenheer, J. A., Booksh, K. (2003). Fluorescence Excitation-Emission Matrix Regional Integration to Quantify Spectra for Dissolved Organic Matter. Environmental Science & Technology 37(24), 5701-5710. \doi{10.1021/es034354c}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @family fluorescence regions
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' regions <- eemsom_regions(som, k = 4, dim = c(8, 8))
#' eemsom_fri(eem_list, regions)
eemsom_fri <- function(eem_list, regions){
  if(!inherits(regions, "eemsom_regions")) stop("regions needs to be of class eemsom_regions!")
  mat <- eem_unfold(eem_list)
  ex <- attr(mat, "ex")
  em <- attr(mat, "em")
  reg <- regions$regions
  same_wl <- isTRUE(all.equal(ex, unique(reg$ex))) && isTRUE(all.equal(em, unique(reg$em)))
  if(!same_wl) stop("The wavelengths of eem_list differ from the regions! Please consider eem_red2smallest() or eem_extend2largest().")
  area <- rep(eem_cellwidth(ex), each = length(em)) * rep(eem_cellwidth(em), times = length(ex))
  rs <- sort(unique(stats::na.omit(reg$region)))
  total_area <- sum(area[!is.na(reg$region)])
  vols <- sapply(rs, function(r){
    idx <- which(reg$region == r)
    x <- mat[, idx, drop = FALSE]
    # area of each wavelength pair, missing values do not count
    a <- matrix(area[idx], nrow = nrow(x), ncol = length(idx), byrow = TRUE)
    a[is.na(x)] <- NA
    rowSums(x * a, na.rm = TRUE) * total_area / rowSums(a, na.rm = TRUE)
  })
  vols <- matrix(vols, nrow = nrow(mat))
  vols[!is.finite(vols)] <- NA
  colnames(vols) <- paste0("region_", rs)
  data.frame(sample = rownames(mat), vols, stringsAsFactors = FALSE, row.names = NULL)
}

#' Width of the cells around wavelengths
#'
#' @param wl numeric vector of wavelengths
#'
#' @return numeric vector, distance between the midpoints to the neighbouring wavelengths
#' @noRd
eem_cellwidth <- function(wl){
  if(length(wl) == 1) return(1)
  d <- diff(wl)
  (c(d[1], d) + c(d, d[length(d)])) / 2
}

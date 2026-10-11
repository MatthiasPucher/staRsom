#' Cluster the units of a self-organising map
#'
#' @description The prototypes (codebook vectors) of the units are grouped by hierarchical clustering or k-means. Samples belong to the cluster of their unit. This two-level approach (Vesanto and Alhoniemi 2000) is faster and more robust than clustering the samples directly.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param k number of clusters, see \code{\link{eemsom_kselect}} for help choosing k
#' @param method \code{"complete"} (default, Ejarque-Gonzalez and Butturini 2014) or any other method of \code{\link[stats]{hclust}}, or \code{"kmeans"}
#' @param seed seed for the random number generator, only used for \code{method = "kmeans"}
#'
#' @return object of class eemsom with the cluster of each unit in the element \code{clusters}
#' @export
#'
#' @importFrom stats hclust cutree kmeans dist setNames
#'
#' @references Vesanto, J., Alhoniemi, E. (2000). Clustering of the self-organizing map. IEEE Transactions on Neural Networks 11(3), 586-600. \doi{10.1109/72.846731}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @seealso \code{\link{eemsom_kselect}}, \code{\link{eemsom_split}}
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' som <- eemsom_cluster(som, k = 3)
#' eemsom_samples(som)
#' eemsom_plot_map(som, type = "cluster")
eemsom_cluster <- function(som, k, method = "complete", seed = 42){
  stopifnot(inherits(som, "eemsom"))
  som$clusters <- eemsom_clusterunits(som$som$codes[[1]], k = k, method = method, seed = seed)
  som$cluster_method <- method
  som
}

#' Help choosing the number of clusters
#'
#' @description The units of a self-organising map are clustered with different numbers of clusters k. Several criteria are calculated for each k.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param k numeric vector of numbers of clusters to be tested
#' @param method clustering method, see \code{\link{eemsom_cluster}}
#' @param seed seed for the random number generator, only used for \code{method = "kmeans"}
#'
#' @details
#' \itemize{
#'   \item silhouette: mean silhouette width of the units (Rousseeuw 1987), higher is better
#'   \item negative: number of units with negative silhouette width, i.e. units that are closer to another cluster, lower is better (Ejarque-Gonzalez and Butturini 2014)
#'   \item davies_bouldin: Davies-Bouldin index, lower is better
#'   \item wss: within-cluster sum of squares, look for an elbow
#' }
#'
#' @return data frame with one row per k
#' @export
#'
#' @importFrom cluster silhouette
#'
#' @references Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. Journal of Computational and Applied Mathematics 20, 53-65. \doi{10.1016/0377-0427(87)90125-7}
#'
#' Davies, D. L., Bouldin, D. W. (1979). A Cluster Separation Measure. IEEE Transactions on Pattern Analysis and Machine Intelligence 1(2), 224-227. \doi{10.1109/TPAMI.1979.4766909}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' eemsom_kselect(som, k = 2:6)
eemsom_kselect <- function(som, k = 2:9, method = "complete", seed = 42){
  stopifnot(inherits(som, "eemsom"))
  eemsom_kcriteria(som$som$codes[[1]], k = k, method = method, seed = seed)
}

#' Plot silhouette widths of clustered units
#'
#' @description The silhouette width of each unit shows how well it fits into its cluster. Values close to 1 indicate a good fit, negative values indicate that the unit is closer to another cluster. Ejarque-Gonzalez and Butturini (2014) chose the number of fluorescence components by the mean silhouette width and the number of negative silhouettes.
#'
#' @param x object of class eemsom with clusters (see \code{\link{eemsom_cluster}}) or of class eemsom_regions (see \code{\link{eemsom_regions}})
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @references Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. Journal of Computational and Applied Mathematics 20, 53-65. \doi{10.1016/0377-0427(87)90125-7}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' som <- eemsom_cluster(som, k = 3)
#' eemsom_plot_silhouette(som)
eemsom_plot_silhouette <- function(x){
  if(inherits(x, "eemsom")){
    if(is.null(x$clusters)) stop("The units have not been clustered yet!")
    codes <- x$som$codes[[1]]
    cl <- x$clusters
  } else if(inherits(x, "eemsom_regions")){
    codes <- x$rsom$codes[[1]]
    cl <- x$unit_regions
  } else {
    stop("x needs to be of class eemsom or eemsom_regions!")
  }
  if(length(unique(cl)) < 2) stop("At least 2 clusters are needed for silhouettes!")
  sil <- silhouette(cl, dist(codes))
  table <- data.frame(unit = seq_along(cl), cluster = factor(sil[, "cluster"]), silhouette = sil[, "sil_width"]) %>%
    dplyr::arrange(cluster, dplyr::desc(silhouette)) %>%
    mutate(pos = seq_along(unit))
  ggplot(table, aes(x = pos, y = silhouette, fill = cluster)) +
    geom_col() +
    geom_hline(yintercept = mean(table$silhouette), linetype = "dashed") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    labs(x = "units", y = "silhouette width")
}

#' Split an eemlist according to the clusters or units of a self-organising map
#'
#' @description The samples are divided into groups that can be analysed separately, e.g. by \code{\link[staRdom]{eem_parafac}}.
#'
#' @param eem_list object of class eemlist, containing samples used to calculate the self-organising map. Usually, these are the original (not normalised) EEMs.
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param by \code{"cluster"} (default) or \code{"unit"}
#' @param min_samples minimum number of samples per group. A warning is shown for smaller groups because a PARAFAC model might not be stable.
#'
#' @details Samples are matched by their names. Groups without samples are omitted.
#'
#' @return a named list of eemlists, names are "cluster_1", "cluster_2", ... or "unit_1", ...
#' @export
#'
#' @seealso \code{\link{eemsom_cluster}}, \code{\link[staRdom]{eem_parafac}}
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' som <- eemsom_cluster(som, k = 3)
#' eem_groups <- eemsom_split(eem_list, som, min_samples = 3)
#' lengths(eem_groups)
#'
#' \donttest{
#' # separate PARAFAC models for each group
#' # (with real data sets, groups need enough samples)
#' pf_groups <- lapply(eem_groups, eem_parafac, comps = 2, nstart = 5, cores = 1)
#' }
eemsom_split <- function(eem_list, som, by = "cluster", min_samples = 20){
  stopifnot(inherits(som, "eemsom"))
  if(!inherits(eem_list, "eemlist")) stop("eem_list is not of class eemlist!")
  by <- match.arg(by, c("cluster", "unit"))
  if(by == "cluster" & is.null(som$clusters)) stop("The units have not been clustered yet!")
  smp <- eemsom_samples(som)
  names <- eem_names(eem_list)
  if(!all(names %in% smp$sample)) stop("The following samples were not used in the self-organising map: ", paste(setdiff(names, smp$sample), collapse = ", "), ". Please use eemsom_predict() for new samples.")
  group <- smp[[by]][match(names, smp$sample)]
  groups <- sort(unique(group))
  res <- lapply(groups, function(g){
    eem_list[group == g] %>%
      `class<-`("eemlist")
  }) %>%
    setNames(paste0(by, "_", groups))
  small <- lengths(res) < min_samples
  if(any(small)) warning("The following groups contain less than ", min_samples, " samples, PARAFAC models might not be stable: ", paste(names(res)[small], collapse = ", "))
  res
}

#' Map new samples to a self-organising map
#'
#' @description New EEMs are normalised like the training data and assigned to their best matching unit.
#'
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param eem_list object of class eemlist with the same excitation and emission wavelengths as the training data
#'
#' @details Missing values in the new samples are ignored when calculating the distances.
#'
#' @return data frame with columns sample, unit, x, y, distance and cluster like \code{\link{eemsom_samples}}
#' @export
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, rlen = 100)
#' som <- eemsom_cluster(som, k = 3)
#' # usually, new samples are mapped. Here, the training samples are used.
#' eemsom_predict(som, eem_list)
eemsom_predict <- function(som, eem_list){
  stopifnot(inherits(som, "eemsom"))
  eempf_som_unsupported(som, "eemsom_predict")
  mat <- eem2mat(eem_list, normalise = "none", na = "keep")
  if(!isTRUE(all.equal(attr(mat, "ex"), som$ex)) | !isTRUE(all.equal(attr(mat, "em"), som$em))) stop("The wavelengths of eem_list differ from the self-organising map! Please consider eem_red2smallest() or eem_extend2largest().")
  mat <- mat[, som$keep, drop = FALSE]
  mat <- mat / eem_normfactors(mat, som$normalise)
  d <- eemsom_distmat(mat, som$som$codes[[1]])
  unit <- apply(d, 1, which.min)
  cluster <- if(is.null(som$clusters)) rep(NA, length(unit)) else som$clusters[unit]
  data.frame(sample = rownames(mat), unit = unit, x = som$som$grid$pts[unit, 1], y = som$som$grid$pts[unit, 2],
             distance = sqrt(d[cbind(seq_along(unit), unit)]), cluster = cluster, stringsAsFactors = FALSE)
}

#' Cluster codebook vectors
#'
#' @param codes matrix of codebook vectors
#' @param k number of clusters
#' @param method hclust method or "kmeans"
#' @param seed seed for kmeans
#'
#' @return integer vector of clusters
#' @noRd
eemsom_clusterunits <- function(codes, k, method = "complete", seed = 42){
  if(!is.numeric(k) | length(k) != 1 || k < 1 || k > nrow(codes)) stop("k needs to be a number between 1 and the number of units (", nrow(codes), ")!")
  hc_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  method <- match.arg(method, c(hc_methods, "kmeans"))
  if(method == "kmeans"){
    cl <- with_seed(seed, kmeans(codes, centers = k, nstart = 25)$cluster)
  } else {
    cl <- cutree(hclust(dist(codes), method = method), k = k)
  }
  unname(cl)
}

#' Criteria for numbers of clusters
#'
#' @param codes matrix of codebook vectors
#' @param k numbers of clusters
#' @param method hclust method or "kmeans"
#' @param seed seed for kmeans
#'
#' @return data frame
#' @noRd
eemsom_kcriteria <- function(codes, k, method = "complete", seed = 42){
  k <- k[k >= 2 & k < nrow(codes)]
  if(length(k) == 0) stop("k needs to be at least 2 and smaller than the number of units (", nrow(codes), ")!")
  d <- dist(codes)
  lapply(k, function(kk){
    cl <- eemsom_clusterunits(codes, k = kk, method = method, seed = seed)
    sil <- silhouette(cl, d)[, "sil_width"]
    data.frame(k = kk, silhouette = mean(sil), negative = sum(sil < 0), davies_bouldin = eemsom_dbindex(codes, cl), wss = eemsom_wss(codes, cl))
  }) %>%
    bind_rows()
}

#' Davies-Bouldin index
#'
#' @param mat matrix, objects in rows
#' @param cl clusters
#'
#' @return numeric
#' @noRd
eemsom_dbindex <- function(mat, cl){
  groups <- sort(unique(cl))
  centres <- t(vapply(groups, function(g) colMeans(mat[cl == g, , drop = FALSE]), numeric(ncol(mat))))
  scatter <- vapply(seq_along(groups), function(i){
    mean(sqrt(rowSums(sweep(mat[cl == groups[i], , drop = FALSE], 2, centres[i, ])^2)))
  }, numeric(1))
  cdist <- as.matrix(dist(centres))
  mean(vapply(seq_along(groups), function(i){
    max((scatter[i] + scatter[-i]) / cdist[i, -i])
  }, numeric(1)))
}

#' Within-cluster sum of squares
#'
#' @param mat matrix, objects in rows
#' @param cl clusters
#'
#' @return numeric
#' @noRd
eemsom_wss <- function(mat, cl){
  sum(vapply(unique(cl), function(g){
    m <- mat[cl == g, , drop = FALSE]
    sum(sweep(m, 2, colMeans(m))^2)
  }, numeric(1)))
}

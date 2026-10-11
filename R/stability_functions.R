#' Stability of a self-organising map
#'
#' @description Self-organising maps are calculated repeatedly from resampled data with the same settings as the original map. The stability of the quantisation and of the neighbourhood relations of the samples is assessed (de Bodt et al. 2002).
#'
#' @param eem_list object of class eemlist, the samples used to calculate \code{som}
#' @param som object of class eemsom, see \code{\link{eemsom}}
#' @param type \code{"bootstrap"} (default, de Bodt et al. 2002) or \code{"loo"} (leave-one-out, Ejarque-Gonzalez and Butturini 2014)
#' @param nrep number of maps (\code{type = "bootstrap"}) or number of maps per left-out sample (\code{type = "loo"})
#' @param radius two samples are neighbours if their units are not more than \code{radius} units apart
#' @param seed seed for the random number generator to get reproducible results
#' @param cores number of parallel calculations
#' @param verbose print infos
#'
#' @details For each map, all samples are mapped to their best matching unit, including samples that were not part of the resampled data.
#'
#' Quantisation: the intra-class sum of squares (SS intra), the sum of the squared distances between samples and their prototypes, is calculated for each map. A low coefficient of variation indicates a stable quantisation.
#'
#' Neighbourhood: for each pair of samples, the stability is the share of maps in which both samples are neighbours. Stable maps have values close to 0 or 1. On a randomly organised map, two samples are neighbours with probability p (element \code{random}), the stability then follows a binomial distribution.
#'
#' With \code{type = "loo"}, \code{nrep} times the number of samples maps are calculated, which can take some time.
#'
#' @return object of class eemsom_stability, a list containing
#' \itemize{
#'   \item quantisation: data frame with columns rep, left_out (only for \code{type = "loo"}) and ss_intra
#'   \item neighbourhood: data frame with columns sample1, sample2 and stability
#'   \item random: probability of two samples being neighbours on a random map
#'   \item nmaps, type, nrep, radius, seed: settings and number of maps
#' }
#' @export
#'
#' @importFrom parallel makeCluster stopCluster parLapply
#'
#' @references de Bodt, E., Cottrell, M., Verleysen, M. (2002). Statistical tools to assess the reliability of self-organizing maps. Neural Networks 15(8-9), 967-978. \doi{10.1016/S0893-6080(02)00071-0}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @seealso \code{\link{eemsom_plot_stability}}
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, dim = c(3, 3), rlen = 20)
#' stab <- eemsom_stability(eem_list, som, nrep = 5)
#' stab
#' eemsom_plot_stability(stab)
#'
#' \donttest{
#' # with real data sets, use more repetitions
#' som <- eemsom(eem_list)
#' stab <- eemsom_stability(eem_list, som, nrep = 20, cores = 2)
#' eemsom_plot_stability(stab)
#' }
eemsom_stability <- function(eem_list, som, type = "bootstrap", nrep = 20, radius = 1, seed = 42, cores = 1, verbose = FALSE){
  stopifnot(inherits(som, "eemsom"))
  eempf_som_unsupported(som, "eemsom_stability")
  type <- match.arg(type, c("bootstrap", "loo"))
  if(is.null(som$settings)) stop("som does not contain its settings, please calculate it again with eemsom()!")
  mat <- eem2mat(eem_list, normalise = som$normalise, na = som$na)
  if(!identical(rownames(mat), som$samples) || !identical(attr(mat, "keep"), som$keep)) stop("eem_list needs to contain the samples used to calculate som, in the same order!")
  n <- nrow(mat)
  tasks <- with_seed(seed, {
    if(type == "bootstrap"){
      lapply(seq_len(nrep), function(r) list(rep = r, left_out = NA_character_, idx = sample.int(n, n, replace = TRUE)))
    } else {
      lapply(seq_len(n * nrep), function(r){
        i <- (r - 1) %/% nrep + 1
        list(rep = r, left_out = rownames(mat)[i], idx = seq_len(n)[-i])
      })
    }
  })
  task_seeds <- with_seed(seed, sample.int(.Machine$integer.max, length(tasks)))
  for(t in seq_along(tasks)) tasks[[t]]$seed <- task_seeds[t]
  settings <- som$settings
  grid <- settings$grid
  if(settings$init == "pca"){
    for(t in seq_along(tasks)) tasks[[t]]$init <- eemsom_pcainit(mat[tasks[[t]]$idx, , drop = FALSE], grid)
  }
  args <- settings[setdiff(names(settings), "init")]
  if(verbose) cat(paste0(length(tasks), " maps are calculated on ", cores, " core", ifelse(cores > 1, "s", ""), "."), fill = TRUE)
  if(cores > 1){
    # detach the worker from the package namespace, so the workers only need kohonen
    worker <- eemsom_stabilityworker
    environment(worker) <- globalenv()
    cl <- makeCluster(min(cores, length(tasks)), type = "PSOCK")
    on.exit(stopCluster(cl))
    codes <- parLapply(cl, tasks, worker, mat = mat, args = args)
  } else {
    codes <- with_seed(seed, lapply(tasks, eemsom_stabilityworker, mat = mat, args = args))
  }
  udist <- unit.distances(grid)
  near <- udist <= radius + 1e-6
  ss_intra <- numeric(length(codes))
  pairs <- utils::combn(n, 2)
  neigh <- numeric(ncol(pairs))
  for(t in seq_along(codes)){
    d <- eemsom_distmat(mat, codes[[t]])
    bmu <- apply(d, 1, which.min)
    ss_intra[t] <- sum(d[cbind(seq_len(n), bmu)])
    neigh <- neigh + near[cbind(bmu[pairs[1, ]], bmu[pairs[2, ]])]
  }
  quantisation <- data.frame(rep = seq_along(codes), ss_intra = ss_intra)
  if(type == "loo") quantisation <- data.frame(rep = seq_along(codes), left_out = vapply(tasks, `[[`, character(1), "left_out"), ss_intra = ss_intra)
  res <- list(quantisation = quantisation,
              neighbourhood = data.frame(sample1 = rownames(mat)[pairs[1, ]], sample2 = rownames(mat)[pairs[2, ]], stability = neigh / length(codes), stringsAsFactors = FALSE),
              random = mean(rowSums(near)) / nrow(udist),
              nmaps = length(codes), type = type, nrep = nrep, radius = radius, seed = seed)
  class(res) <- "eemsom_stability"
  res
}

#' Train one map for the stability assessment
#'
#' @param task list with idx, seed and optionally init
#' @param mat data matrix
#' @param args arguments for kohonen::som
#'
#' @return matrix of codebook vectors
#' @noRd
eemsom_stabilityworker <- function(task, mat, args){
  args$X <- mat[task$idx, , drop = FALSE]
  if(!is.null(task$init)) args$init <- task$init
  set.seed(task$seed)
  do.call(kohonen::som, args)$codes[[1]]
}

#' @rdname eemsom_stability
#' @param x object of class eemsom_stability
#' @param ... further arguments, not used
#' @importFrom stats median sd
#' @export
print.eemsom_stability <- function(x, ...){
  q <- x$quantisation$ss_intra
  st <- x$neighbourhood$stability
  cat("Stability of a self-organising map (", x$type, ", ", x$nmaps, " maps)\n", sep = "")
  cat(" SS intra: mean ", signif(mean(q), 4), ", median ", signif(median(q), 4), ", CV ", signif(sd(q) / mean(q) * 100, 3), " %\n", sep = "")
  cat(" neighbourhood (radius ", x$radius, "): pairs always or never neighbours: ", signif(mean(st == 0 | st == 1) * 100, 3), " %\n", sep = "")
  cat("   expected on a random map: ", signif(((1 - x$random)^x$nmaps + x$random^x$nmaps) * 100, 3), " %\n", sep = "")
  invisible(x)
}

#' Plot the neighbourhood stability of a self-organising map
#'
#' @description A histogram of the neighbourhood stability of all pairs of samples is compared to the distribution expected on a randomly organised map (de Bodt et al. 2002), as in Figure 4B of Ejarque-Gonzalez and Butturini (2014).
#'
#' @param x object of class eemsom_stability, see \code{\link{eemsom_stability}}
#'
#' @return a ggplot object
#' @export
#'
#' @import ggplot2
#' @importFrom stats dbinom
#'
#' @references de Bodt, E., Cottrell, M., Verleysen, M. (2002). Statistical tools to assess the reliability of self-organizing maps. Neural Networks 15(8-9), 967-978. \doi{10.1016/S0893-6080(02)00071-0}
#'
#' Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' @examples
#' data(eem_list)
#'
#' som <- eemsom(eem_list, dim = c(3, 3), rlen = 20)
#' stab <- eemsom_stability(eem_list, som, nrep = 5)
#' eemsom_plot_stability(stab)
eemsom_plot_stability <- function(x){
  if(!inherits(x, "eemsom_stability")) stop("x needs to be of class eemsom_stability!")
  npairs <- nrow(x$neighbourhood)
  k <- 0:x$nmaps
  observed <- data.frame(stability = k / x$nmaps, pairs = tabulate(round(x$neighbourhood$stability * x$nmaps) + 1, nbins = x$nmaps + 1), type = "observed")
  random <- data.frame(stability = k / x$nmaps, pairs = npairs * dbinom(k, x$nmaps, x$random), type = "random map")
  ggplot(observed, aes(x = stability, y = pairs)) +
    geom_col(aes(fill = type), width = 0.8 / x$nmaps) +
    geom_line(data = random, aes(colour = type)) +
    geom_point(data = random, aes(colour = type)) +
    scale_fill_manual(values = "grey60") +
    theme_minimal() +
    labs(x = "neighbourhood stability", y = "pairs of samples", colour = NULL, fill = NULL)
}

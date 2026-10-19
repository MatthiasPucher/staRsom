test_that("eemsom returns a map of class eemsom", {
  expect_s3_class(som_test, "eemsom")
  expect_s3_class(som_test$som, "kohonen")
  expect_equal(som_test$samples, eemR::eem_names(eem_list))
  expect_lte(nrow(som_test$som$codes[[1]]), length(eem_list))
  expect_output(print(som_test), "Self-organising map of 15 EEMs")
  expect_named(som_test$settings, c("grid", "rlen", "alpha", "mode", "cores", "dist.fcts", "init"))
  # the data are always kept, they are needed by the other functions
  expect_false(is.null(eemsom(eem_list, rlen = 10, keep.data = FALSE, seed = 42)$som$data))
})

test_that("results are reproducible with a seed", {
  s1 <- eemsom(eem_list, rlen = 20, seed = 3)
  s2 <- eemsom(eem_list, rlen = 20, seed = 3)
  s3 <- eemsom(eem_list, rlen = 20, seed = 4)
  expect_identical(s1$som$codes, s2$som$codes)
  expect_false(identical(s1$som$codes, s3$som$codes))
  s4 <- eemsom(eem_list, rlen = 20, init = "random", seed = 3)
  s5 <- eemsom(eem_list, rlen = 20, init = "random", seed = 3)
  expect_identical(s4$som$codes, s5$som$codes)
})

test_that("the random number generator state is restored", {
  set.seed(1)
  r1 <- runif(1)
  set.seed(1)
  som <- eemsom(eem_list, rlen = 10, seed = 99)
  expect_equal(runif(1), r1)
})

test_that("without a seed, the random number stream of the session is used", {
  set.seed(5)
  s1 <- eemsom(eem_list, rlen = 20, init = "random")
  r1 <- runif(1)
  set.seed(5)
  s2 <- eemsom(eem_list, rlen = 20, init = "random")
  r2 <- runif(1)
  expect_identical(s1$som$codes, s2$som$codes)
  expect_equal(r1, r2)
  # the stream advances like with any other function using random numbers
  set.seed(5)
  expect_false(isTRUE(all.equal(runif(1), r1)))
})

test_that("seed = NULL is the default of all functions using random numbers", {
  for(f in list(eemsom, eempf_som, eemsom_regions, eemsom_stability, eemsom_cluster, eemsom_kselect)) expect_null(formals(f)$seed)
})

test_that("map size can be set", {
  som <- eemsom(eem_list, dim = c(3, 2), topo = "rectangular", rlen = 20, seed = 42)
  expect_equal(c(som$som$grid$xdim, som$som$grid$ydim), c(3, 2))
  expect_equal(som$som$grid$topo, "rectangular")
  expect_warning(eemsom(eem_list, dim = c(5, 5), rlen = 20, seed = 42), "more units")
  expect_error(eemsom(eem_list, dim = 4, seed = 42), "dim")
  expect_error(eemsom(eem_list, dim = c("a", "b"), seed = 42), "dim")
})

test_that("automatic map size follows the SOM Toolbox", {
  # 5 sqrt(n) units, side ratio from the standard deviations along the first two principal components
  expect_equal(eemsom_dim(list(d = c(1, 1)), 100), c(7, 7))
  expect_equal(eemsom_dim(list(d = c(4, 1)), 100), c(3, 17))
  # extreme ratios (above the square root of the number of units) are replaced by 1
  expect_equal(eemsom_dim(list(d = c(100, 1)), 100), c(7, 7))
  expect_equal(eemsom_dim(list(d = c(1, 1)), 100, topo = "rectangular"), c(7, 7))
  expect_equal(eemsom_dim(list(d = c(1, 1)), 1000), c(12, 13))
  # not more units than samples
  expect_lte(prod(eemsom_dim(list(d = c(2, 1)), 15)), 16)
  expect_equal(c(som_test$som$grid$xdim, som_test$som$grid$ydim), c(3, 5))
  el3 <- structure(eem_list[1:3], class = "eemlist")
  expect_warning(eemsom(el3, rlen = 10, seed = 42), "more units")
})

test_that("other training modes work", {
  expect_s3_class(eemsom(eem_list, mode = "batch", rlen = 20, seed = 42), "eemsom")
  expect_s3_class(eemsom(eem_list, mode = "pbatch", rlen = 20, cores = 1, seed = 42), "eemsom")
})

test_that("each sample is assigned to exactly one unit", {
  smp <- eemsom_samples(som_test)
  expect_equal(nrow(smp), length(eem_list))
  expect_equal(smp$sample, som_test$samples)
  expect_true(all(smp$unit %in% seq_len(nrow(som_test$som$codes[[1]]))))
  expect_true(all(smp$distance >= 0))
  expect_true(all(is.na(smp$cluster)))
})

test_that("number of samples and U-matrix of the unit are given for each sample", {
  smp <- eemsom_samples(som_test)
  expect_named(smp, c("sample", "unit", "x", "y", "distance", "hits", "umatrix", "cluster"))
  expect_equal(smp$hits, as.vector(table(smp$unit)[as.character(smp$unit)]))
  expect_equal(smp$umatrix, eemsom_umatrix(som_test)[smp$unit])
})

test_that("prototypes are returned as eemlist", {
  codes <- eemsom_codes(som_test)
  expect_s3_class(codes, "eemlist")
  expect_length(codes, nrow(som_test$som$codes[[1]]))
  expect_equal(dim(codes[[1]]$x), dim(eem_list[[1]]$x))
  expect_output(print(codes), "V1")
  codes_dn <- eemsom_codes(som_test, denormalise = TRUE)
  u <- som_test$som$unit.classif[1]
  expect_equal(codes_dn[[u]]$x, codes[[u]]$x * mean(som_test$norm_factors[som_test$som$unit.classif == u]))
})

test_that("quality measures are calculated", {
  q <- eemsom_quality(som_test)
  expect_named(q, c("quantisation_error", "topographic_error", "explained_variance"))
  expect_equal(q$quantisation_error, mean(sqrt(som_test$som$distances)))
  expect_true(q$topographic_error >= 0 & q$topographic_error <= 1)
  expect_true(q$explained_variance > 0 & q$explained_variance <= 1)
})

test_that("topographic error counts diagonal neighbours on rectangular maps", {
  som <- eemsom(eem_list, dim = c(3, 3), topo = "rectangular", rlen = 20, seed = 42)
  d <- eemsom_distmat(som$som$data[[1]], som$som$codes[[1]])
  bmu2 <- t(apply(d, 1, function(r) order(r)[1:2]))
  pts <- som$som$grid$pts
  # Chebyshev distance between the two best matching units
  cheb <- pmax(abs(pts[bmu2[, 1], 1] - pts[bmu2[, 2], 1]), abs(pts[bmu2[, 1], 2] - pts[bmu2[, 2], 2]))
  expect_equal(eemsom_quality(som)$topographic_error, mean(cheb > 1))
})

test_that("distances to the codebook vectors are the squared Euclidean distances", {
  mat <- som_test$som$data[[1]]
  codes <- som_test$som$codes[[1]]
  d <- eemsom_distmat(mat, codes)
  expect_equal(dim(d), c(nrow(mat), nrow(codes)))
  expect_equal(d[2, 3], sum((mat[2, ] - codes[3, ])^2))
  expect_null(dimnames(d))
  # missing values: distance from the available variables, scaled up like in kohonen
  m <- mat[1:2, 1:10]
  m[1, 1:4] <- NA
  d <- eemsom_distmat(m, codes[1:2, 1:10])
  expect_equal(d[1, 1], sum((m[1, 5:10] - codes[1, 5:10])^2) * 10 / 6)
  expect_equal(d[2, 2], sum((m[2, ] - codes[2, 1:10])^2))
})

test_that("missing values are handled", {
  expect_lt(sum(som_scatter$keep), length(som_scatter$keep))
  expect_true(all(is.na(eemsom_codes(som_scatter)[[1]]$x[!som_scatter$keep])))
  # na = "keep": pairs missing in all samples are removed, pairs missing in one sample are kept
  all_missing <- Reduce(`&`, lapply(eem_partial, function(eem) is.na(eem$x)))
  expect_equal(sum(som_keep$keep), sum(!all_missing))
  expect_true(anyNA(som_keep$som$data[[1]]))
  codes <- eemsom_codes(som_keep)
  expect_true(all(is.na(codes[[1]]$x[all_missing])))
  expect_false(anyNA(codes[[1]]$x[!all_missing]))
  expect_equal(eemsom(eem_partial, na = "k", rlen = 10, seed = 42)$keep, som_keep$keep)
})

test_that("distances agree between kohonen and staRsom with missing values", {
  smp <- eemsom_samples(som_keep)
  pr <- eemsom_predict(som_keep, eem_partial)
  expect_equal(pr$unit, smp$unit)
  expect_equal(pr$distance, smp$distance)
  expect_equal(eemsom_quality(som_keep)$quantisation_error, mean(smp$distance))
  p <- eemsom_plot_map(som_keep, type = "quality")
  expect_s3_class(p, "ggplot")
})

test_that("automatic map size does not exceed the number of samples", {
  set.seed(1)
  sizes <- sapply(4:60, function(n){
    mat <- matrix(rnorm(n * 10), nrow = n)
    prod(eemsom_dim(eemsom_pca(mat), n)) <= n
  })
  expect_true(all(sizes))
})

test_that("eemsom returns a map of class eemsom", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 50)
  expect_s3_class(som, "eemsom")
  expect_s3_class(som$som, "kohonen")
  expect_equal(som$samples, eemR::eem_names(eem_list))
  expect_lte(nrow(som$som$codes[[1]]), length(eem_list))
  expect_output(print(som), "Self-organising map of 15 EEMs")
})

test_that("results are reproducible with a seed", {
  data(eem_list, package = "staRdom")
  s1 <- eemsom(eem_list, rlen = 50, seed = 3)
  s2 <- eemsom(eem_list, rlen = 50, seed = 3)
  s3 <- eemsom(eem_list, rlen = 50, seed = 4)
  expect_identical(s1$som$codes, s2$som$codes)
  expect_false(identical(s1$som$codes, s3$som$codes))
  s4 <- eemsom(eem_list, rlen = 50, init = "random", seed = 3)
  s5 <- eemsom(eem_list, rlen = 50, init = "random", seed = 3)
  expect_identical(s4$som$codes, s5$som$codes)
})

test_that("the random number generator state is restored", {
  data(eem_list, package = "staRdom")
  set.seed(1)
  r1 <- runif(1)
  set.seed(1)
  som <- eemsom(eem_list, rlen = 20, seed = 99)
  expect_equal(runif(1), r1)
})

test_that("map size can be set", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, dim = c(3, 2), topo = "rectangular", rlen = 20)
  expect_equal(c(som$som$grid$xdim, som$som$grid$ydim), c(3, 2))
  expect_equal(som$som$grid$topo, "rectangular")
  expect_warning(eemsom(eem_list, dim = c(5, 5), rlen = 20), "more units")
  expect_error(eemsom(eem_list, dim = 4), "dim")
})

test_that("other training modes work", {
  data(eem_list, package = "staRdom")
  expect_s3_class(eemsom(eem_list, mode = "batch", rlen = 20), "eemsom")
  expect_s3_class(eemsom(eem_list, mode = "pbatch", rlen = 20, cores = 1), "eemsom")
})

test_that("each sample is assigned to exactly one unit", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 50)
  smp <- eemsom_samples(som)
  expect_equal(nrow(smp), length(eem_list))
  expect_equal(smp$sample, som$samples)
  expect_true(all(smp$unit %in% seq_len(nrow(som$som$codes[[1]]))))
  expect_true(all(smp$distance >= 0))
  expect_true(all(is.na(smp$cluster)))
})

test_that("prototypes are returned as eemlist", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 50)
  codes <- eemsom_codes(som)
  expect_s3_class(codes, "eemlist")
  expect_length(codes, nrow(som$som$codes[[1]]))
  expect_equal(dim(codes[[1]]$x), dim(eem_list[[1]]$x))
  codes_dn <- eemsom_codes(som, denormalise = TRUE)
  u <- som$som$unit.classif[1]
  expect_equal(codes_dn[[u]]$x, codes[[u]]$x * mean(som$norm_factors[som$som$unit.classif == u]))
})

test_that("quality measures are calculated", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 50)
  q <- eemsom_quality(som)
  expect_named(q, c("quantisation_error", "topographic_error", "explained_variance"))
  expect_equal(q$quantisation_error, mean(sqrt(som$som$distances)))
  expect_true(q$topographic_error >= 0 & q$topographic_error <= 1)
  expect_true(q$explained_variance > 0 & q$explained_variance <= 1)
})

test_that("missing values are handled", {
  data(eem_list, package = "staRdom")
  el <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15, 15, 15, 15))
  som_r <- eemsom(el, rlen = 20)
  expect_lt(sum(som_r$keep), length(som_r$keep))
  som_k <- eemsom(el, na = "keep", rlen = 20)
  expect_equal(sum(som_k$keep), length(som_k$keep))
  expect_true(all(is.na(eemsom_codes(som_r)[[1]]$x[!som_r$keep])))
})

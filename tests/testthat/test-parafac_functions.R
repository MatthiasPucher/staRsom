test_that("eempf_som calculates a map of PARAFAC loadings", {
  expect_s3_class(pfsom_test, "eempf_som")
  expect_s3_class(pfsom_test, "eemsom")
  expect_equal(pfsom_test$samples, rownames(pf4[[1]]$A))
  expect_equal(colnames(pfsom_test$som$codes[[1]]), colnames(pf4[[1]]$A))
  expect_equal(unname(rowSums(pfsom_test$som$data[[1]])), rep(1, nrow(pf4[[1]]$A)))
  expect_output(print(pfsom_test), "PARAFAC model with 6 components")
})

test_that("relative loadings are not treated as class memberships by kohonen", {
  pfsom <- eempf_som(pf4[[1]], rlen = 100)
  expect_equal(pfsom$som$dist.fcts, "sumofsquares")
  expect_gt(eemsom_quality(pfsom)$explained_variance, 0.9)
})

test_that("absolute Fmax values can be used", {
  pfsom <- eempf_som(pf4[[1]], normalise = FALSE, rlen = 20)
  pf <- norm2A(eempf_rescaleBC(pf4[[1]]))
  expect_equal(pfsom$som$data[[1]], pf$A, ignore_attr = TRUE)
  expect_output(print(pfsom), "Fmax")
})

test_that("prototype EEMs are reconstructed from PARAFAC components", {
  codes <- eemsom_codes(pfsom_test)
  expect_s3_class(codes, "eemlist")
  pf <- eempf_rescaleBC(pf4[[1]])
  a <- pfsom_test$som$codes[[1]][2, ]
  expect_equal(codes[[2]]$x, pf$B %*% diag(a) %*% t(pf$C), ignore_attr = TRUE)
  expect_equal(codes[[2]]$ex, as.numeric(rownames(pf$C)))
  expect_s3_class(ggeem(pfsom_test, units = 1:2), "ggplot")
  codes_dn <- eemsom_codes(pfsom_test, denormalise = TRUE)
  u <- pfsom_test$som$unit.classif[1]
  expect_equal(codes_dn[[u]]$x, codes[[u]]$x * mean(pfsom_test$norm_factors[pfsom_test$som$unit.classif == u]))
})

test_that("clustering, splitting and plots work with maps of PARAFAC loadings", {
  pfsom <- eemsom_cluster(pfsom_test, k = 3)
  expect_equal(nrow(eemsom_samples(pfsom)), nrow(pf4[[1]]$A))
  groups <- eemsom_split(eem_list, pfsom, min_samples = 1)
  expect_equal(sum(lengths(groups)), length(eem_list))
  expect_s3_class(eemsom_plot_map(pfsom, type = "cluster"), "ggplot")
  p <- eemsom_plot_cplane(pfsom, comps = c(1, 3))
  expect_equal(levels(p$data$pair), c("Comp.1", "Comp.3"))
  expect_equal(nlevels(eemsom_plot_cplane(pfsom)$data$pair), 6)
  expect_error(eemsom_plot_cplane(pfsom, ex = 250, em = 400), "comps")
  expect_error(eemsom_plot_cplane(pfsom, comps = 7), "comps")
})

test_that("new samples are mapped to maps of PARAFAC loadings", {
  # the samples of eem_list are part of the PARAFAC model, so their loadings are reproduced
  pr <- eemsom_predict(pfsom_test, eem_prepared)
  expect_named(pr, c("sample", "unit", "x", "y", "distance", "hits", "umatrix", "cluster"))
  expect_equal(pr$sample, eemR::eem_names(eem_prepared))
  smp <- eemsom_samples(pfsom_test)
  smp <- smp[match(pr$sample, smp$sample), ]
  expect_lt(max(abs(pr$distance - smp$distance)), 0.01)
  expect_gte(mean(pr$unit == smp$unit), 0.8)
  # arguments for A_missing() can be given, the output of A_missing() is suppressed
  expect_silent(pr2 <- eemsom_predict(pfsom_test, eem_prepared, cores = 1, maxit = 100))
  expect_equal(pr2$unit, pr$unit)
  # absolute loadings
  pfsom_abs <- eempf_som(pf4[[1]], normalise = FALSE, rlen = 20)
  pr_abs <- eemsom_predict(pfsom_abs, eem_prepared)
  smp_abs <- eemsom_samples(pfsom_abs)
  expect_lt(max(abs(pr_abs$distance - smp_abs$distance[match(pr_abs$sample, smp_abs$sample)])), 0.01)
  # wavelengths need to fit the model
  expect_error(eemsom_predict(pfsom_test, eem_list), "wavelengths")
  pfsom_old <- pfsom_test
  pfsom_old$pfmodel <- NULL
  expect_error(eemsom_predict(pfsom_old, eem_prepared), "PARAFAC model")
})

test_that("unsupported functions give an error", {
  expect_error(eemsom_regions(pfsom_test), "not available")
  expect_error(eemsom_stability(eem_list, pfsom_test), "not available")
  expect_error(eempf_som(eem_list), "parafac")
})

test_that("component planes of EEM maps reject comps", {
  expect_error(eemsom_plot_cplane(som_test, comps = 1), "PARAFAC")
})

test_that("the Euclidean distance is used for sum-normalised EEMs", {
  el <- eem_list
  for(n in seq_along(el)) el[[n]]$x <- pmax(el[[n]]$x, 0)
  som <- eemsom(el, normalise = "sum", rlen = 20)
  expect_equal(som$som$dist.fcts, "sumofsquares")
})

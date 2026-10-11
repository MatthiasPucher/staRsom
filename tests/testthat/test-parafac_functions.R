test_that("eempf_som calculates a map of PARAFAC loadings", {
  data(pf_models, package = "staRdom")
  pfsom <- eempf_som(pf4[[1]], rlen = 50)
  expect_s3_class(pfsom, "eempf_som")
  expect_s3_class(pfsom, "eemsom")
  expect_equal(pfsom$samples, rownames(pf4[[1]]$A))
  expect_equal(colnames(pfsom$som$codes[[1]]), colnames(pf4[[1]]$A))
  expect_equal(unname(rowSums(pfsom$som$data[[1]])), rep(1, nrow(pf4[[1]]$A)))
  expect_output(print(pfsom), "PARAFAC model with 6 components")
})

test_that("relative loadings are not treated as class memberships by kohonen", {
  data(pf_models, package = "staRdom")
  pfsom <- eempf_som(pf4[[1]], rlen = 100)
  expect_equal(pfsom$som$dist.fcts, "sumofsquares")
  expect_gt(eemsom_quality(pfsom)$explained_variance, 0.9)
})

test_that("absolute Fmax values can be used", {
  data(pf_models, package = "staRdom")
  pfsom <- eempf_som(pf4[[1]], normalise = FALSE, rlen = 20)
  pf <- norm2A(eempf_rescaleBC(pf4[[1]]))
  expect_equal(pfsom$som$data[[1]], pf$A, ignore_attr = TRUE)
  expect_output(print(pfsom), "Fmax")
})

test_that("prototype EEMs are reconstructed from PARAFAC components", {
  data(pf_models, package = "staRdom")
  pfsom <- eempf_som(pf4[[1]], rlen = 20)
  codes <- eemsom_codes(pfsom)
  expect_s3_class(codes, "eemlist")
  pf <- eempf_rescaleBC(pf4[[1]])
  a <- pfsom$som$codes[[1]][2, ]
  expect_equal(codes[[2]]$x, pf$B %*% diag(a) %*% t(pf$C), ignore_attr = TRUE)
  expect_equal(codes[[2]]$ex, as.numeric(rownames(pf$C)))
  expect_s3_class(ggeem(pfsom, units = 1:2), "ggplot")
})

test_that("clustering, splitting and plots work with maps of PARAFAC loadings", {
  data(pf_models, package = "staRdom")
  data(eem_list, package = "staRdom")
  pfsom <- eemsom_cluster(eempf_som(pf4[[1]], rlen = 20), k = 3)
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

test_that("unsupported functions give an error", {
  data(pf_models, package = "staRdom")
  data(eem_list, package = "staRdom")
  pfsom <- eempf_som(pf4[[1]], rlen = 20)
  expect_error(eemsom_regions(pfsom), "not available")
  expect_error(eemsom_predict(pfsom, eem_list), "not available")
  expect_error(eemsom_stability(eem_list, pfsom), "not available")
  expect_error(eempf_som(eem_list), "parafac")
})

test_that("component planes of EEM maps reject comps", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 20)
  expect_error(eemsom_plot_cplane(som, comps = 1), "PARAFAC")
})

test_that("the Euclidean distance is used for sum-normalised EEMs", {
  data(eem_list, package = "staRdom")
  el <- eem_list
  for(n in seq_along(el)) el[[n]]$x <- pmax(el[[n]]$x, 0)
  som <- eemsom(el, normalise = "sum", rlen = 20)
  expect_equal(som$som$dist.fcts, "sumofsquares")
})

test_that("ggeem plots prototypes of all or selected units", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 20)
  p <- ggeem(som)
  expect_s3_class(p, "ggplot")
  expect_equal(nlevels(p$data$sample), nrow(som$som$codes[[1]]))
  expect_true(all(grepl("\\(n = [0-9]+\\)$", levels(p$data$sample))))
  p2 <- ggeem(som, units = c(2, 1))
  expect_equal(nlevels(p2$data$sample), 2)
  expect_error(ggeem(som, units = 100), "units")
})

test_that("map plots are created for all types", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 20)
  for(type in c("counts", "umatrix", "quality", "samples")){
    expect_s3_class(eemsom_plot_map(som, type = type), "ggplot")
  }
  expect_error(eemsom_plot_map(som, type = "cluster"), "clustered")
  expect_s3_class(eemsom_plot_map(som, type = "samples", colour = seq_along(eem_list), labels = TRUE), "ggplot")
  expect_error(eemsom_plot_map(som, type = "samples", colour = 1:2), "colour")
  p <- eemsom_plot_map(som)
  expect_equal(sum(p$data$value) / 6, length(eem_list))
})

test_that("rectangular maps are plotted", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, dim = c(3, 3), topo = "rectangular", rlen = 20)
  p <- eemsom_plot_map(som, type = "umatrix")
  expect_equal(nrow(p$data), 9 * 4)
})

test_that("u-matrix contains mean distances to neighbours", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, dim = c(2, 2), topo = "rectangular", rlen = 20)
  um <- eemsom_umatrix(som)
  codes <- som$som$codes[[1]]
  expect_length(um, 4)
  # on rectangular maps, kohonen counts diagonal units as neighbours
  expect_equal(um[1], mean(sqrt(c(sum((codes[1, ] - codes[2, ])^2), sum((codes[1, ] - codes[3, ])^2), sum((codes[1, ] - codes[4, ])^2)))))
})

test_that("component planes use the closest wavelengths", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 20)
  p <- eemsom_plot_cplane(som, ex = c(276, 350), em = c(311, 450))
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$pair), c("ex 275 nm, em 310 nm", "ex 350 nm, em 450 nm"))
  expect_error(eemsom_plot_cplane(som, ex = 275, em = c(310, 450)), "same length")
  el <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15, 15, 15, 15))
  som_na <- eemsom(el, rlen = 20)
  expect_error(eemsom_plot_cplane(som_na, ex = 350, em = 350), "removed")
})

test_that("training progress is plotted", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 20)
  p <- eemsom_plot_training(som)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 20)
})

test_that("component planes can be scaled individually", {
  data(eem_list, package = "staRdom")
  som <- eemsom(eem_list, rlen = 20)
  p_common <- eemsom_plot_cplane(som, ex = c(275, 350), em = c(310, 450))
  p_each <- eemsom_plot_cplane(som, ex = c(275, 350), em = c(310, 450), scale = "each")
  ranges <- tapply(p_each$data$value, p_each$data$pair, range)
  expect_true(all(vapply(ranges, function(r) isTRUE(all.equal(r, c(0, 1))), logical(1))))
  expect_equal(p_common$data$value[p_common$data$pair == levels(p_common$data$pair)[1]][1], unname(som$som$codes[[1]][1, "275|310"]))
  expect_error(eemsom_plot_cplane(som, ex = 275, em = 310, scale = "foo"))
  data(pf_models, package = "staRdom")
  pfsom <- eempf_som(pf4[[1]], rlen = 20)
  expect_equal(range(eemsom_plot_cplane(pfsom, comps = 2, scale = "each")$data$value), c(0, 1))
})

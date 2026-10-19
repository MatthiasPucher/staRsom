test_that("ggeem plots prototypes of all or selected units", {
  p <- ggeem(som_test)
  expect_s3_class(p, "ggplot")
  expect_equal(nlevels(p$data$sample), nrow(som_test$som$codes[[1]]))
  expect_true(all(grepl("\\(n = [0-9]+\\)$", levels(p$data$sample))))
  p2 <- ggeem(som_test, units = c(2, 1))
  expect_equal(nlevels(p2$data$sample), 2)
  expect_error(ggeem(som_test, units = 100), "units")
})

test_that("map plots are created for all types", {
  for(type in c("counts", "umatrix", "quality", "samples")){
    expect_s3_class(eemsom_plot_map(som_test, type = type), "ggplot")
  }
  expect_error(eemsom_plot_map(som_test, type = "cluster"), "clustered")
  expect_s3_class(eemsom_plot_map(som_clustered, type = "cluster"), "ggplot")
  expect_s3_class(eemsom_plot_map(som_test, type = "samples", colour = seq_along(eem_list), labels = TRUE), "ggplot")
  expect_error(eemsom_plot_map(som_test, type = "samples", colour = 1:2), "colour")
  p <- eemsom_plot_map(som_test)
  expect_equal(sum(p$data$value) / 6, length(eem_list))
  expect_false(any(sapply(p$layers, function(l) inherits(l$geom, "GeomText"))))
  pu <- eemsom_plot_map(som_test, unit_labels = TRUE)
  text_layers <- Filter(function(l) inherits(l$geom, "GeomText"), pu$layers)
  expect_length(text_layers, 1)
  expect_equal(text_layers[[1]]$data$unit, seq_len(nrow(som_test$som$codes[[1]])))
  expect_s3_class(eemsom_plot_map(som_clustered, type = "samples", labels = TRUE, unit_labels = TRUE), "ggplot")
})

test_that("rectangular maps are plotted", {
  som <- eemsom(eem_list, dim = c(3, 3), topo = "rectangular", rlen = 20, seed = 42)
  p <- eemsom_plot_map(som, type = "umatrix")
  expect_equal(nrow(p$data), 9 * 4)
})

test_that("u-matrix contains mean distances to neighbours", {
  som <- eemsom(eem_list, dim = c(2, 2), topo = "rectangular", rlen = 20, seed = 42)
  um <- eemsom_umatrix(som)
  codes <- som$som$codes[[1]]
  expect_length(um, 4)
  # on rectangular maps, kohonen counts diagonal units as neighbours
  expect_equal(um[1], mean(sqrt(c(sum((codes[1, ] - codes[2, ])^2), sum((codes[1, ] - codes[3, ])^2), sum((codes[1, ] - codes[4, ])^2)))))
})

test_that("component planes use the closest wavelengths", {
  p <- eemsom_plot_cplane(som_test, ex = c(276, 350), em = c(311, 450))
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$pair), c("ex 275 nm, em 310 nm", "ex 350 nm, em 450 nm"))
  expect_error(eemsom_plot_cplane(som_test, ex = 275, em = c(310, 450)), "same length")
  expect_error(eemsom_plot_cplane(som_scatter, ex = 350, em = 350), "removed")
})

test_that("training progress is plotted", {
  p <- eemsom_plot_training(som_small)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 20)
})

test_that("component planes can be scaled individually", {
  p_common <- eemsom_plot_cplane(som_test, ex = c(275, 350), em = c(310, 450))
  p_each <- eemsom_plot_cplane(som_test, ex = c(275, 350), em = c(310, 450), scale = "each")
  ranges <- tapply(p_each$data$value, p_each$data$pair, range)
  expect_true(all(vapply(ranges, function(r) isTRUE(all.equal(r, c(0, 1))), logical(1))))
  expect_equal(p_common$data$value[p_common$data$pair == levels(p_common$data$pair)[1]][1], unname(som_test$som$codes[[1]][1, "275|310"]))
  expect_error(eemsom_plot_cplane(som_test, ex = 275, em = 310, scale = "foo"))
  expect_equal(range(eemsom_plot_cplane(pfsom_test, comps = 2, scale = "each")$data$value), c(0, 1))
})

test_that("discrete scales use viridis", {
  viridis <- function(n) ggplot2::scale_fill_viridis_d()$palette(n)
  fills <- function(p) unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
  expect_setequal(fills(eemsom_plot_map(som_clustered, type = "cluster")), viridis(3))
  expect_setequal(fills(eemsom_plot_silhouette(som_clustered)), viridis(3))
  p <- eemsom_plot_map(som_test, type = "samples", colour = rep(c("a", "b", "c"), 5))
  expect_setequal(unique(ggplot2::ggplot_build(p)$data[[2]]$colour), ggplot2::scale_colour_viridis_d(end = 0.9)$palette(3))
})


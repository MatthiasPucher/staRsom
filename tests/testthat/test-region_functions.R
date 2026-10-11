test_that("regions are found for all wavelength pairs", {
  reg <- eemsom_regions(som_test, dim = c(6, 6), rlen = 20)
  expect_s3_class(reg, "eemsom_regions")
  expect_equal(nrow(reg$regions), length(som_test$keep))
  expect_named(reg$regions, c("ex", "em", "region"))
  expect_equal(reg$kselect$k, 2:9)
  expect_equal(reg$k, reg$kselect$k[which.max(reg$kselect$silhouette)])
  expect_setequal(unique(stats::na.omit(reg$regions$region)), seq_len(reg$k))
  expect_length(reg$unit_regions, 36)
  expect_output(print(reg), "regions:")
})

test_that("number of regions can be set and results are reproducible", {
  r1 <- eemsom_regions(som_test, k = 3, dim = c(5, 5), rlen = 20, seed = 1)
  r2 <- eemsom_regions(som_test, k = 3, dim = c(5, 5), rlen = 20, seed = 1)
  expect_equal(r1$k, 3)
  expect_null(r1$kselect)
  expect_identical(r1$regions, r2$regions)
  expect_error(eemsom_regions(som_test, dim = 5), "dim")
  expect_error(eemsom_regions(som_test, k = 50, dim = c(3, 3), rlen = 5), "R-mode")
  expect_error(eemsom_regions(som_test, k = "a", dim = c(3, 3), rlen = 5), "candidates")
})

test_that("candidates for the number of regions can be given", {
  r <- eemsom_regions(som_test, k = 2:4, dim = c(5, 5), rlen = 20)
  expect_equal(r$kselect$k, 2:4)
  expect_true(r$k %in% 2:4)
  expect_equal(r$k, r$kselect$k[which.max(r$kselect$silhouette)])
  r_null <- eemsom_regions(som_test, k = NULL, dim = c(5, 5), rlen = 20)
  r_default <- eemsom_regions(som_test, dim = c(5, 5), rlen = 20)
  expect_equal(r_null$kselect, r_default$kselect)
  expect_equal(r_default$kselect$k, 2:9)
})

test_that("regions of planes and correlation matrix are similar", {
  el <- eem_exclude(eem_list, list(ex = eem_list[[1]]$ex[c(FALSE, TRUE)], em = eem_list[[1]]$em[c(FALSE, TRUE, TRUE)]))
  som <- eemsom(el, rlen = 50)
  rp <- eemsom_regions(som, k = 2, method = "planes", dim = c(5, 5), rlen = 20)
  rc <- eemsom_regions(som, k = 2, method = "cormat", dim = c(5, 5), rlen = 20)
  tab <- table(rp$regions$region, rc$regions$region)
  # most wavelength pairs are in corresponding regions
  expect_gt(max(sum(diag(tab)), sum(tab) - sum(diag(tab))) / sum(tab), 0.7)
})

test_that("removed wavelength pairs get no region", {
  reg <- eemsom_regions(som_scatter, k = 2, dim = c(5, 5), rlen = 20)
  expect_true(all(is.na(reg$regions$region[!som_scatter$keep])))
  # with na = "keep", pairs missing in all samples get no region either
  reg_keep <- eemsom_regions(som_keep, k = 2, dim = c(5, 5), rlen = 20)
  expect_true(all(is.na(reg_keep$regions$region[!som_keep$keep])))
  expect_setequal(unique(stats::na.omit(reg_keep$regions$region)), 1:2)
})

test_that("regions are plotted", {
  reg <- eemsom_regions(som_test, k = 3, dim = c(5, 5), rlen = 20)
  p <- ggeem(reg)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$data$region, "factor")
  ps <- eemsom_plot_silhouette(reg)
  expect_s3_class(ps, "ggplot")
  expect_equal(nrow(ps$data), 25)
})

test_that("fluorescence regional integration is area normalised", {
  reg <- eemsom_regions(som_test, k = 3, dim = c(5, 5), rlen = 20)
  fri <- eemsom_fri(eem_list, reg)
  expect_equal(fri$sample, eemR::eem_names(eem_list))
  expect_equal(names(fri), c("sample", paste0("region_", 1:3)))
  # a constant EEM gives the same normalised volume in every region
  el <- eem_list
  for(n in seq_along(el)) el[[n]]$x[] <- 2
  fri_c <- eemsom_fri(el, reg)
  vol <- unname(as.matrix(fri_c[, -1]))
  expect_equal(vol[1, ], rep(vol[1, 1], 3))
  # equals intensity x total area of all pairs with a region (regular grid: 5 nm ex x 4 nm em)
  expect_equal(vol[1, 1], 2 * 5 * 4 * sum(!is.na(reg$regions$region)))
  # missing values reduce the area of the region, the normalised volume stays the same
  el[[1]]$x[1:30, 1:10] <- NA
  fri_na <- eemsom_fri(el, reg)
  expect_equal(unname(as.matrix(fri_na[, -1]))[1, ], vol[1, ])
  expect_error(eemsom_fri(eem_range(eem_list, ex = c(250, Inf)), reg), "wavelengths")
  expect_error(eemsom_fri(eem_list, som_test), "eemsom_regions")
})

test_that("cell widths follow the wavelength spacing", {
  expect_equal(eem_cellwidth(c(1, 2, 3)), c(1, 1, 1))
  expect_equal(eem_cellwidth(c(0, 2, 6)), c(2, 3, 4))
  expect_equal(eem_cellwidth(5), 1)
})

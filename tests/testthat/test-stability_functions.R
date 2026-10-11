test_that("bootstrap stability is calculated", {
  st <- eemsom_stability(eem_list, som_small, nrep = 3)
  expect_s3_class(st, "eemsom_stability")
  expect_equal(st$nmaps, 3)
  expect_equal(nrow(st$quantisation), 3)
  expect_true(all(st$quantisation$ss_intra > 0))
  n <- length(eem_list)
  expect_equal(nrow(st$neighbourhood), n * (n - 1) / 2)
  expect_true(all(st$neighbourhood$stability %in% c(0, 1/3, 2/3, 1)))
  expect_output(print(st), "SS intra")
  expect_output(eemsom_stability(eem_list, som_small, nrep = 2, verbose = TRUE), "maps calculated")
})

test_that("random neighbourhood probability fits the grid", {
  som <- eemsom(eem_list, dim = c(2, 2), topo = "rectangular", rlen = 20)
  st <- eemsom_stability(eem_list, som, nrep = 2)
  # all units of a 2 x 2 rectangular map are neighbours
  expect_equal(st$random, 1)
  expect_true(all(st$neighbourhood$stability == 1))
})

test_that("stability is reproducible and independent of the number of cores", {
  skip_on_cran()
  s1 <- eemsom_stability(eem_list, som_small, nrep = 3, seed = 2)
  s2 <- eemsom_stability(eem_list, som_small, nrep = 3, seed = 2)
  expect_identical(s1, s2)
  s3 <- eemsom_stability(eem_list, som_small, nrep = 3, seed = 2, cores = 2)
  expect_equal(s3$quantisation, s1$quantisation)
  expect_equal(s3$neighbourhood, s1$neighbourhood)
})

test_that("settings of the map are used", {
  som <- eemsom(eem_list, dim = c(3, 3), rlen = 10, init = "random", radius = 1)
  expect_equal(som$settings$radius, 1)
  st <- eemsom_stability(eem_list, som, nrep = 2)
  expect_equal(st$nmaps, 2)
  st_keep <- eemsom_stability(eem_partial, som_keep, nrep = 2)
  expect_equal(st_keep$nmaps, 2)
})

test_that("leave-one-out stability is calculated", {
  st <- eemsom_stability(eem_list, som_small, type = "loo", nrep = 1)
  expect_equal(st$nmaps, length(eem_list))
  expect_equal(st$quantisation$left_out, eemR::eem_names(eem_list))
})

test_that("wrong input is rejected", {
  el <- eem_list[c(2, 1, 3:15)]
  class(el) <- "eemlist"
  expect_error(eemsom_stability(el, som_small, nrep = 2), "same order")
  som <- som_small
  som$settings <- NULL
  expect_error(eemsom_stability(eem_list, som, nrep = 2), "settings")
})

test_that("stability is plotted", {
  st <- eemsom_stability(eem_list, som_small, nrep = 3)
  p <- eemsom_plot_stability(st)
  expect_s3_class(p, "ggplot")
  expect_equal(sum(p$data$pairs), nrow(st$neighbourhood))
  expect_error(eemsom_plot_stability(som_small), "eemsom_stability")
})

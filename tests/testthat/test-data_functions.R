test_that("eem2mat unfolds an eemlist", {
  data(eem_list, package = "staRdom")
  m <- eem2mat(eem_list)
  expect_equal(dim(m), c(length(eem_list), length(eem_list[[1]]$ex) * length(eem_list[[1]]$em)))
  expect_equal(rownames(m), eemR::eem_names(eem_list))
  expect_equal(colnames(m)[1], paste0(eem_list[[1]]$ex[1], "|", eem_list[[1]]$em[1]))
  expect_equal(unname(apply(m, 1, sd)), rep(1, length(eem_list)))
  expect_length(attr(m, "norm_factors"), length(eem_list))
})

test_that("normalisation methods work", {
  data(eem_list, package = "staRdom")
  expect_equal(unname(apply(eem2mat(eem_list, normalise = "max"), 1, max)), rep(1, length(eem_list)))
  expect_equal(unname(rowSums(eem2mat(eem_list, normalise = "sum"))), rep(1, length(eem_list)))
  expect_equal(attr(eem2mat(eem_list, normalise = "none"), "norm_factors"), rep(1, length(eem_list)))
  expect_error(eem2mat(eem_list, normalise = "foo"))
})

test_that("mat2eem reverses eem2mat", {
  data(eem_list, package = "staRdom")
  m <- eem2mat(eem_list)
  back <- mat2eem(m, norm_factors = attr(m, "norm_factors"))
  expect_s3_class(back, "eemlist")
  expect_equal(eemR::eem_names(back), eemR::eem_names(eem_list))
  for(n in seq_along(eem_list)) expect_equal(back[[n]]$x, eem_list[[n]]$x, ignore_attr = TRUE)
  expect_equal(back[[1]]$ex, eem_list[[1]]$ex)
  expect_equal(back[[1]]$em, eem_list[[1]]$em)
})

test_that("missing values are removed and restored as NA", {
  data(eem_list, package = "staRdom")
  el <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15, 15, 15, 15))
  m <- eem2mat(el)
  expect_false(anyNA(m))
  expect_lt(ncol(m), length(attr(m, "keep")))
  expect_true(anyNA(eem2mat(el, na = "keep")))
  back <- mat2eem(m, norm_factors = attr(m, "norm_factors"))
  x0 <- el[[1]]$x
  x1 <- back[[1]]$x
  expect_true(all(is.na(x1[is.na(x0)])))
  expect_equal(x1[!is.na(x1)], x0[!is.na(x1)])
})

test_that("wrong input is rejected", {
  data(eem_list, package = "staRdom")
  expect_error(eem2mat(unclass(eem_list)), "eemlist")
  el <- eem_list
  el[[2]]$ex <- el[[2]]$ex + 1
  expect_error(eem2mat(el), "dimensions mismatch")
  m <- eem2mat(eem_list)
  expect_error(mat2eem(m[, -1], ex = attr(m, "ex"), em = attr(m, "em"), keep = attr(m, "keep")), "columns")
})

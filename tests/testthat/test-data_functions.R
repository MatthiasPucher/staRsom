test_that("eem2mat unfolds an eemlist", {
  m <- eem2mat(eem_list)
  expect_equal(dim(m), c(length(eem_list), length(eem_list[[1]]$ex) * length(eem_list[[1]]$em)))
  expect_equal(rownames(m), eemR::eem_names(eem_list))
  expect_equal(colnames(m)[1], paste0(eem_list[[1]]$ex[1], "|", eem_list[[1]]$em[1]))
  expect_equal(unname(apply(m, 1, sd)), rep(1, length(eem_list)))
  expect_named(attr(m, "norm_factors"), eemR::eem_names(eem_list))
})

test_that("normalisation methods work", {
  expect_equal(unname(apply(eem2mat(eem_list, normalise = "max"), 1, max)), rep(1, length(eem_list)))
  expect_equal(unname(rowSums(eem2mat(eem_list, normalise = "sum"))), rep(1, length(eem_list)))
  expect_equal(unname(attr(eem2mat(eem_list, normalise = "none"), "norm_factors")), rep(1, length(eem_list)))
  expect_error(eem2mat(eem_list, normalise = "foo"))
})

test_that("normalisation factors need to be positive", {
  el <- eem_list
  el[[1]]$x <- -abs(el[[1]]$x) - 1
  expect_error(eem2mat(el, normalise = "max"), "not positive")
  expect_error(eem2mat(el, normalise = "sum"), "not positive")
  el[[1]]$x[] <- 0
  expect_error(eem2mat(el, normalise = "sd"), "not positive")
})

test_that("mat2eem reverses eem2mat and gives printable EEMs", {
  m <- eem2mat(eem_list)
  back <- mat2eem(m, norm_factors = attr(m, "norm_factors"))
  expect_s3_class(back, "eemlist")
  expect_equal(eemR::eem_names(back), eemR::eem_names(eem_list))
  for(n in seq_along(eem_list)) expect_equal(back[[n]]$x, eem_list[[n]]$x, ignore_attr = TRUE)
  expect_equal(back[[1]]$ex, eem_list[[1]]$ex)
  expect_equal(back[[1]]$em, eem_list[[1]]$em)
  expect_output(print(back), eemR::eem_names(eem_list)[1])
  expect_output(print(summary(back)), "is_blank_corrected")
  expect_length(mat2eem(m[1, ], ex = attr(m, "ex"), em = attr(m, "em"), keep = attr(m, "keep")), 1)
})

test_that("missing values are removed and restored as NA", {
  m <- eem2mat(eem_scatter)
  expect_false(anyNA(m))
  expect_lt(ncol(m), length(attr(m, "keep")))
  back <- mat2eem(m, norm_factors = attr(m, "norm_factors"))
  x0 <- eem_scatter[[1]]$x
  x1 <- back[[1]]$x
  expect_true(all(is.na(x1[is.na(x0)])))
  expect_equal(x1[!is.na(x1)], x0[!is.na(x1)])
})

test_that("na = 'keep' keeps pairs that are missing in some samples only", {
  full <- do.call(rbind, lapply(eem_partial, function(eem) as.vector(eem$x)))
  m <- eem2mat(eem_partial, na = "keep")
  expect_true(anyNA(m))
  expect_true(all(colSums(!is.na(m)) > 0))
  expect_equal(attr(m, "keep"), colSums(!is.na(full)) > 0)
  expect_equal(sum(!attr(m, "keep")), sum(colSums(is.na(full)) == nrow(full)))
  expect_equal(attr(eem2mat(eem_partial, na = "k"), "keep"), attr(m, "keep"))
})

test_that("wrong input is rejected", {
  expect_error(eem2mat(unclass(eem_list)), "eemlist")
  el <- eem_list
  el[[2]]$ex <- el[[2]]$ex + 1
  expect_error(eem2mat(el), "dimensions mismatch")
  el <- eem_list
  el[[3]]$x <- t(el[[3]]$x)
  expect_error(eem2mat(el), "dimensions mismatch")
  el <- eem_list
  el[[2]]$sample <- el[[1]]$sample
  expect_error(eem2mat(el), "unique")
  m <- eem2mat(eem_list)
  expect_error(mat2eem(m[, -1], ex = attr(m, "ex"), em = attr(m, "em"), keep = attr(m, "keep")), "columns")
})

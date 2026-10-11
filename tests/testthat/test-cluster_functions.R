test_that("units are clustered", {
  expect_s3_class(som_clustered, "eemsom")
  expect_length(som_clustered$clusters, nrow(som_test$som$codes[[1]]))
  expect_setequal(unique(som_clustered$clusters), 1:3)
  expect_equal(eemsom_samples(som_clustered)$cluster, som_clustered$clusters[som_test$som$unit.classif])
  expect_output(print(som_clustered), "clusters: 3")
  expect_equal(eemsom_cluster(som_test, k = 3, method = "kmeans", seed = 1)$clusters, eemsom_cluster(som_test, k = 3, method = "kmeans", seed = 1)$clusters)
  expect_setequal(unique(eemsom_cluster(som_test, k = 2, method = "ward.D2")$clusters), 1:2)
  expect_error(eemsom_cluster(som_test, k = 100), "k needs")
  expect_error(eemsom_cluster(som_test, k = 2, method = "foo"))
})

test_that("criteria for the number of clusters are calculated", {
  ks <- eemsom_kselect(som_test, k = 2:5)
  expect_named(ks, c("k", "silhouette", "negative", "davies_bouldin", "wss"))
  expect_equal(ks$k, 2:5)
  expect_true(all(ks$silhouette >= -1 & ks$silhouette <= 1))
  expect_true(all(diff(ks$wss) <= 0))
  expect_equal(eemsom_kselect(som_test, k = c(1, 3, 100))$k, 3)
  expect_error(eemsom_kselect(som_test, k = 100), "k needs")
})

test_that("Davies-Bouldin index is correct for a simple case", {
  mat <- matrix(c(0, 0, 1, 0, 10, 0, 11, 0), ncol = 2, byrow = TRUE)
  # scatter 0.5 in both clusters, distance of centres 10
  expect_equal(eemsom_dbindex(mat, c(1, 1, 2, 2)), 0.1)
  expect_equal(eemsom_wss(mat, c(1, 1, 2, 2)), 1)
  # identical centres are not separated
  expect_equal(eemsom_dbindex(matrix(c(0, 0, 0, 0, 5, 5), ncol = 2, byrow = TRUE), c(1, 2, 3)), Inf)
})

test_that("silhouettes are plotted", {
  expect_error(eemsom_plot_silhouette(som_test), "clustered")
  p <- eemsom_plot_silhouette(som_clustered)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(som_test$som$codes[[1]]))
  expect_error(eemsom_plot_silhouette(eem_list), "class")
})

test_that("eemlist is split into a complete partition", {
  groups <- suppressWarnings(eemsom_split(eem_list, som_clustered))
  expect_type(groups, "list")
  expect_true(all(vapply(groups, inherits, logical(1), "eemlist")))
  expect_equal(sort(unlist(lapply(groups, eemR::eem_names), use.names = FALSE)), sort(eemR::eem_names(eem_list)))
  expect_equal(sum(lengths(groups)), length(eem_list))
  expect_true(all(grepl("^cluster_", names(groups))))
  expect_warning(eemsom_split(eem_list, som_clustered, min_samples = 20), "less than 20")
  expect_no_warning(eemsom_split(eem_list, som_clustered, min_samples = 1))
  units <- eemsom_split(eem_list, som_clustered, by = "unit", min_samples = 1)
  expect_equal(length(units), length(unique(som_test$som$unit.classif)))
  expect_error(eemsom_split(eem_list, som_test), "clustered")
  expect_error(eemsom_split(unclass(eem_list), som_clustered), "eemlist")
})

test_that("split rejects samples not used in the map", {
  el <- eem_list
  el[[1]]$sample <- "new_sample"
  expect_error(eemsom_split(el, som_clustered, min_samples = 1), "new_sample")
})

test_that("new samples are mapped like training samples", {
  pr <- eemsom_predict(som_clustered, eem_list)
  smp <- eemsom_samples(som_clustered)
  expect_equal(pr$unit, smp$unit)
  expect_equal(pr$distance, smp$distance)
  expect_equal(pr$cluster, smp$cluster)
  el <- eem_list
  el[[1]]$x <- el[[1]]$x * 5
  expect_equal(eemsom_predict(som_clustered, el)$unit[1], smp$unit[1])
  el2 <- eem_range(eem_list, ex = c(250, Inf))
  expect_error(eemsom_predict(som_clustered, el2), "wavelengths")
  # missing values in new samples are ignored, the distance is scaled up
  el3 <- eem_list
  el3[[1]]$x[1:5, 1:3] <- NA
  pr3 <- eemsom_predict(som_clustered, el3)
  expect_equal(pr3$unit[-1], smp$unit[-1])
  expect_false(is.na(pr3$distance[1]))
})

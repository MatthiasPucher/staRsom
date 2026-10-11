# fixtures shared by the tests, calculated once
data(eem_list, package = "staRdom")
data(pf_models, package = "staRdom")

# EEMs with scatter removed (missing in all samples) and additionally missing values in one sample
eem_scatter <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15, 15, 15, 15))
eem_partial <- eem_scatter
eem_partial[[1]]$x[20:40, 10:20] <- NA

som_test <- eemsom(eem_list, rlen = 50)
som_clustered <- eemsom_cluster(som_test, k = 3)
som_small <- eemsom(eem_list, dim = c(3, 3), rlen = 20)
som_scatter <- eemsom(eem_scatter, rlen = 20)
som_keep <- eemsom(eem_partial, na = "keep", rlen = 20)
pfsom_test <- eempf_som(pf4[[1]], rlen = 20)

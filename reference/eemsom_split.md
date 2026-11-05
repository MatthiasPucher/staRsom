# Split an eemlist according to the clusters or units of a self-organising map

The samples are divided into groups that can be analysed separately,
e.g. by
[`eem_parafac`](https://rdrr.io/pkg/staRdom/man/eem_parafac.html).

## Usage

``` r
eemsom_split(eem_list, som, by = "cluster", min_samples = 20)
```

## Arguments

- eem_list:

  object of class eemlist, containing samples used to calculate the
  self-organising map. Usually, these are the original (not normalised)
  EEMs.

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

- by:

  `"cluster"` (default) or `"unit"`

- min_samples:

  minimum number of samples per group. A warning is shown for smaller
  groups because a PARAFAC model might not be stable.

## Value

a named list of eemlists, names are "cluster_1", "cluster_2", ... or
"unit_1", ...

## Details

Samples are matched by their names. Groups without samples are omitted.

## See also

[`eemsom_cluster`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md),
[`eem_parafac`](https://rdrr.io/pkg/staRdom/man/eem_parafac.html)

Other clustering:
[`eemsom_cluster()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md),
[`eemsom_kselect()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_kselect.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md)

Other map results:
[`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md),
[`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md),
[`eemsom_quality()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md),
[`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
som <- eemsom_cluster(som, k = 3)
eem_groups <- eemsom_split(eem_list, som, min_samples = 3)
#> Warning: The following groups contain less than 3 samples, PARAFAC models might not be stable: cluster_3
lengths(eem_groups)
#> cluster_1 cluster_2 cluster_3 
#>         6         8         1 

# \donttest{
# separate PARAFAC models for the groups with at least 5 samples
# (with real data sets, groups need more samples, starts and iterations)
eem_groups <- eem_groups[lengths(eem_groups) >= 5]
pf_groups <- lapply(eem_groups, eem_parafac, comps = 2, nstart = 1, maxit = 500, cores = 1)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> Warning: Calculating the 2 components model, 0 out of 1 models converged! You might want to increase the number of initialisations (nstart) or iterations (maxit).
#> Warning: The PARAFAC model with 2 components did not converge! Increasing the number of initialisations (nstart) or iterations (maxit) might solve the problem.
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
# }
```

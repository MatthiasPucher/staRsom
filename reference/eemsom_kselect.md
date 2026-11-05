# Help choosing the number of clusters

The units of a self-organising map are clustered with different numbers
of clusters k. Several criteria are calculated for each k.

## Usage

``` r
eemsom_kselect(som, k = 2:9, method = "complete", seed = NULL)
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

- k:

  numeric vector of numbers of clusters to be tested

- method:

  clustering method, see
  [`eemsom_cluster`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md)

- seed:

  seed for the random number generator, only used for
  `method = "kmeans"`. With `NULL` (default), the current random number
  stream of the session is used, so
  [`set.seed`](https://rdrr.io/r/base/Random.html) before the call gives
  reproducible results. With a number, the seed is set for this
  calculation only and the random number state of the session is
  restored afterwards.

## Value

data frame with one row per k

## Details

- silhouette: mean silhouette width of the units (Rousseeuw 1987),
  higher is better

- negative: number of units with negative silhouette width, i.e. units
  that are closer to another cluster, lower is better (Ejarque-Gonzalez
  and Butturini 2014)

- davies_bouldin: Davies-Bouldin index, lower is better (Inf if two
  clusters have identical centres)

- wss: within-cluster sum of squares, look for an elbow

## References

Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the
interpretation and validation of cluster analysis. Journal of
Computational and Applied Mathematics 20, 53-65.
[doi:10.1016/0377-0427(87)90125-7](https://doi.org/10.1016/0377-0427%2887%2990125-7)

Davies, D. L., Bouldin, D. W. (1979). A Cluster Separation Measure. IEEE
Transactions on Pattern Analysis and Machine Intelligence 1(2), 224-227.
[doi:10.1109/TPAMI.1979.4766909](https://doi.org/10.1109/TPAMI.1979.4766909)

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

Other clustering:
[`eemsom_cluster()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
eemsom_kselect(som, k = 2:6)
#>   k silhouette negative davies_bouldin       wss
#> 1 2  0.5203208        0      0.6639666 8333.2481
#> 2 3  0.5702871        0      0.5061423 2676.1187
#> 3 4  0.5125841        0      0.3559383 1863.7990
#> 4 5  0.3985742        0      0.5463554 1220.6489
#> 5 6  0.4120996        0      0.4334395  863.8047
```

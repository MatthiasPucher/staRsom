# Assignment of samples to units of a self-organising map

Assignment of samples to units of a self-organising map

## Usage

``` r
eemsom_samples(som)
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

## Value

data frame with columns sample, unit (number of the best matching unit),
x and y (position of the unit on the map), distance (Euclidean distance
between sample and prototype), hits (number of samples mapped to the
unit), umatrix (mean distance between the prototype of the unit and the
prototypes of its neighbours, see
[`eemsom_plot_map`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md))
and cluster (NA if no clustering was done)

## Details

With missing values (`na = "keep"`), the distances are calculated from
the available wavelength pairs of a sample and scaled up to the number
of all wavelength pairs, as in kohonen.

The columns hits and umatrix help to find unusual samples (outliers):
such samples often occupy a unit alone (hits = 1) that is separated from
its neighbours by a high U-matrix value. Their distance to the prototype
can be small though, because the prototype adapts to them. The
leave-one-out distance calculated by
[`eemsom_stability`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
with `type = "loo"` is a more reliable indicator.

## See also

Other map results:
[`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md),
[`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md),
[`eemsom_quality()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md),
[`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
eemsom_samples(som)
#>          sample unit   x         y  distance hits  umatrix cluster
#> 1      d0680sfK    9 3.5 2.5980762  2.300073    1 36.48733      NA
#> 2       d1266sf    3 3.5 0.8660254  9.054573    2 10.92376      NA
#> 3      d1268sfK    3 3.5 0.8660254 11.834547    2 10.92376      NA
#> 4      d1543sfK   13 1.5 4.3301270  7.938452    5 11.58416      NA
#> 5      d1549sfK   13 1.5 4.3301270  9.879895    5 11.58416      NA
#> 6      d1550sfK    7 1.5 2.5980762  4.606082    3 15.87562      NA
#> 7      d1557sfK   13 1.5 4.3301270  4.190952    5 11.58416      NA
#> 8      d1559sfK   13 1.5 4.3301270  8.586989    5 11.58416      NA
#> 9      d1565sfK    7 1.5 2.5980762  6.557108    3 15.87562      NA
#> 10     d1567sfK   13 1.5 4.3301270 18.078294    5 11.58416      NA
#> 11     d1568sfK    7 1.5 2.5980762  8.106521    3 15.87562      NA
#> 12       d492sf    1 1.5 0.8660254  3.276772    2 21.19327      NA
#> 13       d661sf    1 1.5 0.8660254  3.092985    2 21.19327      NA
#> 14       d667sf    2 2.5 0.8660254  1.694963    1 14.83993      NA
#> 15 dSFB445pwd1t   15 3.5 4.3301270  2.842423    1 55.31539      NA
```

# Quality measures of a self-organising map

Quality measures of a self-organising map

## Usage

``` r
eemsom_quality(som)
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

## Value

data frame with one row

## Details

- quantisation_error: mean Euclidean distance between samples and their
  prototypes, see
  [`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)
  for the handling of missing values

- topographic_error: share of samples whose best and second best
  matching units are not neighbours (six neighbours on hexagonal maps,
  eight including the diagonals on rectangular maps, see
  [`unit.distances`](https://rdrr.io/pkg/kohonen/man/unit.distances.html))

- explained_variance: share of the total sum of squares of the
  (normalised) data explained by the prototypes

## References

Kohonen, T. (2013). Essentials of the self-organizing map. Neural
Networks 37, 52-65.
[doi:10.1016/j.neunet.2012.09.018](https://doi.org/10.1016/j.neunet.2012.09.018)

## See also

Other map results:
[`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md),
[`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md),
[`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md),
[`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
eemsom_quality(som)
#>   quantisation_error topographic_error explained_variance
#> 1           6.802709                 0          0.9473782
```

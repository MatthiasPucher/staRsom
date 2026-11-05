# Map new samples to a self-organising map

New EEMs are assigned to their best matching unit. For maps of EEMs
([`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)),
they are normalised like the training data. For maps of PARAFAC loadings
([`eempf_som`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)),
their loadings are calculated first with the fixed emission and
excitation loadings of the PARAFAC model
([`A_missing`](https://rdrr.io/pkg/staRdom/man/A_missing.html)).

## Usage

``` r
eemsom_predict(som, eem_list, cores = 1, ...)
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)
  or
  [`eempf_som`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)

- eem_list:

  object of class eemlist with the same excitation and emission
  wavelengths as the training data or the PARAFAC model

- cores:

  number of cores used by
  [`A_missing`](https://rdrr.io/pkg/staRdom/man/A_missing.html) (maps of
  PARAFAC loadings only)

- ...:

  additional arguments passed on to
  [`A_missing`](https://rdrr.io/pkg/staRdom/man/A_missing.html) and
  [`eem_parafac`](https://rdrr.io/pkg/staRdom/man/eem_parafac.html)
  (maps of PARAFAC loadings only), e.g. `maxit`. By default,
  `nstart = 1` is used because the loadings with fixed B and C modes do
  not depend on the start.

## Value

data frame with columns sample, unit, x, y, distance, hits, umatrix and
cluster like
[`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md).
hits is the number of training samples mapped to the unit. A distance
much higher than the distances of the training samples indicates a
sample that is not represented by the map.

## Details

Maps of EEMs: only the wavelength pairs used in the map are considered.
Missing values in the new samples are ignored when calculating the
distances, which are scaled up to the number of all wavelength pairs as
in kohonen (see
[`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)).

Maps of PARAFAC loadings: the new EEMs need the same preparation as the
EEMs of the PARAFAC model (e.g. scatter removal). Their Fmax loadings
are scaled like in
[`eempf_som`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)
(relative or absolute) and the distances refer to these loadings.

## See also

Other map results:
[`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md),
[`eemsom_quality()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md),
[`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md),
[`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
som <- eemsom_cluster(som, k = 3)
# usually, new samples are mapped. Here, the training samples are used.
eemsom_predict(som, eem_list)
#>          sample unit   x         y  distance hits  umatrix cluster
#> 1      d0680sfK    9 3.5 2.5980762  2.300073    1 36.48733       1
#> 2       d1266sf    3 3.5 0.8660254  9.054573    2 10.92376       1
#> 3      d1268sfK    3 3.5 0.8660254 11.834547    2 10.92376       1
#> 4      d1543sfK   13 1.5 4.3301270  7.938452    5 11.58416       2
#> 5      d1549sfK   13 1.5 4.3301270  9.879895    5 11.58416       2
#> 6      d1550sfK    7 1.5 2.5980762  4.606082    3 15.87562       2
#> 7      d1557sfK   13 1.5 4.3301270  4.190952    5 11.58416       2
#> 8      d1559sfK   13 1.5 4.3301270  8.586989    5 11.58416       2
#> 9      d1565sfK    7 1.5 2.5980762  6.557108    3 15.87562       2
#> 10     d1567sfK   13 1.5 4.3301270 18.078294    5 11.58416       2
#> 11     d1568sfK    7 1.5 2.5980762  8.106521    3 15.87562       2
#> 12       d492sf    1 1.5 0.8660254  3.276772    2 21.19327       1
#> 13       d661sf    1 1.5 0.8660254  3.092985    2 21.19327       1
#> 14       d667sf    2 2.5 0.8660254  1.694963    1 14.83993       1
#> 15 dSFB445pwd1t   15 3.5 4.3301270  2.842423    1 55.31539       3

# map of PARAFAC loadings: the loadings of the new samples are calculated first
data(pf_models)
data(eem_list_outliers)
pfsom <- eempf_som(pf4[[1]], rlen = 100, seed = 42)
new_samples <- eem_list_outliers |>
  eem_rem_scat(remove_scatter = c(TRUE, TRUE, TRUE, TRUE),
               remove_scatter_width = c(15, 15, 18, 19)) |>
  eem_interp(type = 1, extend = FALSE, cores = 1)
eemsom_predict(pfsom, new_samples)
#>       sample unit   x        y   distance hits    umatrix cluster
#> 1 dsfb676psp   46 1.0 8.660254 0.06066783    7 0.02492000      NA
#> 2  dsgb447wt   25 5.5 4.330127 0.05285394    6 0.01610523      NA
```

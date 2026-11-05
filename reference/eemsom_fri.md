# Fluorescence regional integration

The fluorescence of each sample is integrated within each region.
Following Chen et al. (2003), the volumes are normalised by the area of
the region, so regions of different size can be compared.
Ejarque-Gonzalez and Butturini (2014) quantified their fluorescence
components from the self-organising map in this way.

## Usage

``` r
eemsom_fri(eem_list, regions)
```

## Arguments

- eem_list:

  object of class eemlist with the same wavelengths as the
  self-organising map used for the regions. Usually, these are the
  original (not normalised) EEMs.

- regions:

  object of class eemsom_regions, see
  [`eemsom_regions`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)

## Value

data frame with column sample and one column per region (region_1,
region_2, ...) containing the area-normalised volumes

## Details

The volume of region i is \\\Phi_i = \sum I(ex, em) \Delta ex \Delta
em\\ over all wavelength pairs of the region. The area-normalised volume
is \\\Phi_i \cdot A / A_i\\, with A being the total area of all regions
and \\A_i\\ the area of region i. Missing values in a sample are ignored
and the area of the region is reduced accordingly.

## References

Chen, W., Westerhoff, P., Leenheer, J. A., Booksh, K. (2003).
Fluorescence Excitation-Emission Matrix Regional Integration to Quantify
Spectra for Dissolved Organic Matter. Environmental Science & Technology
37(24), 5701-5710.
[doi:10.1021/es034354c](https://doi.org/10.1021/es034354c)

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

Other fluorescence regions:
[`eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md),
[`ggeem.eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
regions <- eemsom_regions(som, k = 4, dim = c(8, 8), seed = 42)
eemsom_fri(eem_list, regions)
#>          sample    region_1 region_2  region_3   region_4
#> 1      d0680sfK -946.199953 7300.221  7632.905 10941.5058
#> 2       d1266sf  872.158976 8038.719 10946.700  9047.7166
#> 3      d1268sfK    7.757121 7933.168  4679.323 10145.2257
#> 4      d1543sfK    9.558215 3061.654 11349.520 12488.3607
#> 5      d1549sfK    7.421150 2363.786  8167.564  8910.0850
#> 6      d1550sfK    5.112928 2455.799  3095.346 11301.4860
#> 7      d1557sfK    7.443168 2403.499  9622.169  9966.6183
#> 8      d1559sfK    7.314715 1597.210  6590.773  7078.7004
#> 9      d1565sfK    3.156784 1795.550  1061.570  8776.8975
#> 10     d1567sfK    8.411932 2947.562 17018.794 18144.9301
#> 11     d1568sfK    4.848648 2287.318  3942.471 13017.9332
#> 12       d492sf  429.832122 5970.442  8948.428  8874.6858
#> 13       d661sf  471.578490 5313.502  8076.655  7926.7696
#> 14       d667sf  572.856939 6704.365  8923.545  8649.1114
#> 15 dSFB445pwd1t 6684.784958 1945.030  5452.151   426.5932
```

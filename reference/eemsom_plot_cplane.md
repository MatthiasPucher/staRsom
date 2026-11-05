# Plot component planes of a self-organising map

A component plane shows the value of one variable in the prototypes of
all units. Similar component planes indicate correlated variables. For
maps of EEMs, the variables are wavelength pairs, for maps of PARAFAC
loadings
([`eempf_som`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)),
the variables are the PARAFAC components.

## Usage

``` r
eemsom_plot_cplane(som, ex = NULL, em = NULL, comps = NULL, scale = "common")
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)
  or
  [`eempf_som`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)

- ex:

  numeric vector of excitation wavelengths (maps of EEMs only)

- em:

  numeric vector of emission wavelengths, same length as `ex`. Each pair
  of `ex` and `em` gives one component plane.

- comps:

  numeric vector or names of PARAFAC components (maps of PARAFAC
  loadings only). NULL (default) plots all components.

- scale:

  `"common"` (default) uses one colour scale for all component planes,
  so values can be compared between planes. `"each"` scales each
  component plane from 0 (minimum) to 1 (maximum), which shows patterns
  of weak peaks or components more clearly.

## Value

a ggplot object

## Details

The closest wavelengths available in the data are used. With
`scale = "common"`, values are on the scale of the normalised data.

Ejarque-Gonzalez and Butturini (2014) compared the component planes of
all wavelength pairs by their correlation to discriminate fluorescence
components, see
[`eemsom_regions`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md).

## References

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

Other plots:
[`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md),
[`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md),
[`ggeem.eemsom()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md),
[`ggeem.eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
# fluorescence peaks B, T, A and C (Coble 1996)
eemsom_plot_cplane(som, ex = c(275, 275, 260, 350), em = c(310, 340, 450, 450))

eemsom_plot_cplane(som, ex = c(275, 275, 260, 350), em = c(310, 340, 450, 450), scale = "each")


data(pf_models)
pfsom <- eempf_som(pf4[[1]], rlen = 100, seed = 42)
eemsom_plot_cplane(pfsom, scale = "each")
```

# Plot prototype EEMs of a self-organising map

The prototype EEMs (codebook vectors) of all units or selected units are
plotted. If all units are plotted, the facets are arranged like the
units on the map (top row of the map first).

## Usage

``` r
# S3 method for class 'eemsom'
ggeem(data, fill_max = FALSE, units = NULL, denormalise = FALSE, ...)
```

## Arguments

- data:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

- fill_max:

  sets the maximum fluorescence value for the colour scale, see
  [`ggeem`](https://rdrr.io/pkg/staRdom/man/ggeem.html)

- units:

  optional numeric vector of units to be plotted

- denormalise:

  logical, see
  [`eemsom_codes`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md)

- ...:

  parameters passed on to
  [`ggeem`](https://rdrr.io/pkg/staRdom/man/ggeem.html), e.g. `colpal`
  or `contour`

## Value

a ggplot object

## Details

The facet titles show the unit and the number of samples mapped to this
unit. On a hexagonal map, every second row is shifted by half a unit,
which cannot be shown by the facets.

Ejarque-Gonzalez and Butturini (2014) used the prototype EEMs to
describe the characteristic fluorescence of the regions of the map.

## References

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

[`ggeem`](https://rdrr.io/pkg/staRdom/man/ggeem.html),
[`eemsom_codes`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md)

Other plots:
[`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md),
[`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md),
[`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md),
[`ggeem.eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
ggeem(som)

ggeem(som, units = c(1, 2), contour = TRUE)
```

# Plot a self-organising map

Different properties of the units are shown on the map.

## Usage

``` r
eemsom_plot_map(
  som,
  type = "counts",
  colour = NULL,
  labels = FALSE,
  unit_labels = FALSE
)
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

- type:

  property to be shown:

  - `"counts"`: number of samples mapped to each unit

  - `"umatrix"`: unified distance matrix, mean Euclidean distance
    between the prototype of a unit and the prototypes of its direct
    neighbours. High values separate groups of similar units.

  - `"quality"`: mean distance between the samples of a unit and its
    prototype, see
    [`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)

  - `"cluster"`: clusters of the units

  - `"samples"`: samples as points on their units, the units are
    coloured by cluster if available

- colour:

  optional vector with one value per sample (e.g. sampling site from
  metadata) used to colour the samples with `type = "samples"`

- labels:

  logical, whether sample names are shown with `type = "samples"`

- unit_labels:

  logical, whether the numbers of the units are shown. They are used by
  [`ggeem.eemsom`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md),
  [`eemsom_codes`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md)
  and
  [`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md).

## Value

a ggplot object

## Details

Neighbours are defined like in
[`unit.distances`](https://rdrr.io/pkg/kohonen/man/unit.distances.html):
six units on hexagonal maps, eight units including the diagonals on
rectangular maps.

With `type = "samples"`, samples of the same unit are arranged in a
circle around the centre of the unit.

Ejarque-Gonzalez and Butturini (2014) used the U-matrix to find regions
of similar units and projected the samples, coloured by sampling site
and hydrological conditions, onto the map to explore the patterns among
samples.

## References

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

Other plots:
[`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md),
[`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md),
[`ggeem.eemsom()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md),
[`ggeem.eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
eemsom_plot_map(som, unit_labels = TRUE)

eemsom_plot_map(som, type = "umatrix")

eemsom_plot_map(som, type = "samples", colour = substr(eem_names(eem_list), 1, 2), labels = TRUE)
```

# Plot the training progress of a self-organising map

The mean distance between samples and their closest unit during the
training is plotted. The curve should flatten at the end of the
training, otherwise `rlen` might be increased.

## Usage

``` r
eemsom_plot_training(som)
```

## Arguments

- som:

  object of class eemsom, see
  [`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)

## Value

a ggplot object

## See also

Other plots:
[`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md),
[`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md),
[`ggeem.eemsom()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md),
[`ggeem.eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
eemsom_plot_training(som)
```

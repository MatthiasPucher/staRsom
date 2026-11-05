# Plot fluorescence regions

The regions found by
[`eemsom_regions`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)
are shown in the excitation-emission space.

## Usage

``` r
# S3 method for class 'eemsom_regions'
ggeem(data, fill_max = FALSE, ...)
```

## Arguments

- data:

  object of class eemsom_regions

- fill_max:

  not used, for consistency with
  [`ggeem`](https://rdrr.io/pkg/staRdom/man/ggeem.html)

- ...:

  further arguments, not used

## Value

a ggplot object

## See also

Other fluorescence regions:
[`eemsom_fri()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_fri.md),
[`eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)

Other plots:
[`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md),
[`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md),
[`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md),
[`ggeem.eemsom()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
regions <- eemsom_regions(som, dim = c(8, 8), seed = 42)
ggeem(regions)
```

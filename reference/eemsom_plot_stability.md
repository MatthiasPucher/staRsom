# Plot the neighbourhood stability of a self-organising map

A histogram of the neighbourhood stability of all pairs of samples is
compared to the distribution expected on a randomly organised map (de
Bodt et al. 2002), as in Figure 4B of Ejarque-Gonzalez and Butturini
(2014).

## Usage

``` r
eemsom_plot_stability(x)
```

## Arguments

- x:

  object of class eemsom_stability, see
  [`eemsom_stability`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)

## Value

a ggplot object

## References

de Bodt, E., Cottrell, M., Verleysen, M. (2002). Statistical tools to
assess the reliability of self-organizing maps. Neural Networks 15(8-9),
967-978.
[doi:10.1016/S0893-6080(02)00071-0](https://doi.org/10.1016/S0893-6080%2802%2900071-0)

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

Other stability:
[`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)

Other plots:
[`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md),
[`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
[`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md),
[`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md),
[`ggeem.eemsom()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md),
[`ggeem.eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, dim = c(3, 3), rlen = 20, seed = 42)
stab <- eemsom_stability(eem_list, som, nrep = 3, seed = 42)
eemsom_plot_stability(stab)
```

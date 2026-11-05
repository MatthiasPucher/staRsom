# Package index

## Self-organising maps

Train a map from EEMs and look at its samples, prototypes and quality.

- [`eemsom()`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)
  [`print(`*`<eemsom>`*`)`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)
  : Calculate a self-organising map from EEM data
- [`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)
  : Assignment of samples to units of a self-organising map
- [`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md)
  : Prototype EEMs of a self-organising map
- [`eemsom_quality()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md)
  : Quality measures of a self-organising map
- [`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md)
  : Map new samples to a self-organising map

## Plots

- [`ggeem(`*`<eemsom>`*`)`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom.md)
  : Plot prototype EEMs of a self-organising map
- [`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md)
  : Plot a self-organising map
- [`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md)
  : Plot component planes of a self-organising map
- [`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md)
  : Plot the training progress of a self-organising map

## Clusters and groups of samples

- [`eemsom_kselect()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_kselect.md)
  : Help choosing the number of clusters
- [`eemsom_cluster()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md)
  : Cluster the units of a self-organising map
- [`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md)
  : Plot silhouette widths of clustered units
- [`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)
  : Split an eemlist according to the clusters or units of a
  self-organising map

## Fluorescence regions

- [`eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)
  [`print(`*`<eemsom_regions>`*`)`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)
  : Find fluorescence regions from correlated component planes
- [`ggeem(`*`<eemsom_regions>`*`)`](https://matthiaspucher.github.io/staRsom/reference/ggeem.eemsom_regions.md)
  : Plot fluorescence regions
- [`eemsom_fri()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_fri.md)
  : Fluorescence regional integration

## Stability and outliers

- [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  [`print(`*`<eemsom_stability>`*`)`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  : Stability of a self-organising map
- [`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md)
  : Plot the neighbourhood stability of a self-organising map

## PARAFAC models

- [`eempf_som()`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)
  [`print(`*`<eempf_som>`*`)`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)
  : Calculate a self-organising map from the loadings of a PARAFAC model

## Data conversion

- [`eem2mat()`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md)
  : Unfold EEMs into a matrix of samples and wavelength pairs
- [`mat2eem()`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md)
  : Fold a matrix of unfolded EEMs back into an eemlist

## Package

- [`staRsom`](https://matthiaspucher.github.io/staRsom/reference/staRsom-package.md)
  [`staRsom-package`](https://matthiaspucher.github.io/staRsom/reference/staRsom-package.md)
  : staRsom: Self-Organising Maps for EEMs from DOM

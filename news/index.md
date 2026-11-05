# Changelog

## staRsom 0.3.0

- Prepared for CRAN.
- No fixed seed in functions: `seed = NULL` is the default of
  [`eemsom()`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md),
  [`eempf_som()`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md),
  [`eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md),
  [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md),
  [`eemsom_cluster()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md)
  and
  [`eemsom_kselect()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_kselect.md).
  The random number stream of the session is used, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  gives reproducible results. A number sets the seed for this
  calculation only, `seed = 42` gives the results of earlier versions.
- [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  with `cores = 1` keeps the random number stream of the session;
  before, the seed of its last map remained set.
- No CITATION file: `citation("staRsom")` is generated from DESCRIPTION.
  For the package from CRAN, it contains the DOI of the package. The
  method paper of Ejarque-Gonzalez and Butturini (2014) is cited in
  DESCRIPTION, README, help pages and vignette.
- DESCRIPTION links to the GitHub repository, the README shows the logo,
  the installation from GitHub and how the package was developed.
- Shorter examples of
  [`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)
  and
  [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md).
- Website built with pkgdown:
  <https://matthiaspucher.github.io/staRsom/>.

## staRsom 0.1.3

- The vignette shows the staRsom logo, derived from the staRdom logo.

## staRsom 0.1.2

- [`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)
  and
  [`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md)
  give the number of samples on the unit (`hits`) and its U-matrix value
  (`umatrix`) for each sample. Unusual samples often occupy a unit alone
  that is separated from its neighbours.
- [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  with `type = "loo"` returns the leave-one-out distance of each sample
  (element `loo`): the distance to its best matching unit on the maps
  calculated without it. High values indicate candidates for outliers in
  a subsequent PARAFAC analysis;
  [`print()`](https://rdrr.io/r/base/print.html) shows the highest
  values.
- All plots use viridis colour scales, also for discrete values:
  clusters and sample colours in
  [`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
  regions in [`ggeem()`](https://rdrr.io/pkg/staRdom/man/ggeem.html),
  silhouettes and the random map in
  [`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md).
  Numeric sample colours in `eemsom_plot_map(type = "samples")`
  (e.g. leave-one-out distances) use the continuous viridis scale.
  Sample points omit the light yellow end of viridis to remain visible
  on the white background.
- Automatic map size: rounding could lead to more units than samples
  (e.g. 18 units for 17 samples), as documented the number of units is
  now limited to the number of samples.
- Vignette: new section on outliers. Outliers are usually kept in the
  map; its indicators are compared with the leverage of a PARAFAC model.
- Vignette rewritten in the style of the staRdom vignettes: practical
  advice addressed to the reader, task-oriented sections, an explanation
  of self-organising maps for ecologists and biogeochemists, a diagram
  of the analysis steps and notes on parallel processing.

## staRsom 0.1.1

- [`eem2mat()`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md)
  always removes wavelength pairs that are missing in all samples, also
  with `na = "keep"`. kohonen filled them with random values that were
  never trained, so prototype EEMs, component planes and regions were
  wrong there.
- EEMs created by
  [`mat2eem()`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md)
  (e.g. by
  [`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md))
  carry the attributes of eemR and can be printed.
- Distances of samples with missing values are scaled like in kohonen,
  so
  [`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md),
  [`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md),
  [`eemsom_quality()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md)
  and
  [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  agree. Distances are calculated by matrix products, which is much
  faster for large data sets.
- [`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md)
  gains `unit_labels` to show the numbers of the units.
- [`eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md):
  `k` can be a vector of candidates (default `2:9`) or a single number.
- [`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md)
  maps new samples to maps of PARAFAC loadings
  ([`eempf_som()`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)):
  their loadings are calculated with the fixed emission and excitation
  loadings of the model
  ([`A_missing()`](https://rdrr.io/pkg/staRdom/man/A_missing.html) of
  staRdom). `eempf_som` objects store the PARAFAC model in the element
  `pfmodel` instead of `components`.
- Automatic map size: extreme side ratios are replaced by 1 as in the
  SOM Toolbox.
- Input checks: normalisation factors need to be positive, sample names
  need to be unique, abbreviations of `na` are accepted, `keep.data` is
  always TRUE.
- [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  shows the progress with `verbose = TRUE`. Its results differ slightly
  from version 0.1.0 because the seeds of the repetitions are drawn
  differently.
- Documentation of run times and advice for larger data sets in
  [`eemsom()`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md),
  [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  and the vignette; cross-references between related functions;
  `Language: en-GB`.

## staRsom 0.1.0

First release.

- [`eem2mat()`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md)
  and
  [`mat2eem()`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md)
  convert eemlists to matrices and back, including normalisation and
  handling of missing values.
- [`eemsom()`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)
  calculates a self-organising map from EEM data;
  [`eemsom_samples()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md),
  [`eemsom_codes()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md)
  and
  [`eemsom_quality()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md)
  extract results.
- [`ggeem()`](https://rdrr.io/pkg/staRdom/man/ggeem.html) method for
  `eemsom` objects plots the prototype EEMs arranged like the map;
  [`eemsom_plot_map()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_map.md),
  [`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md)
  and
  [`eemsom_plot_training()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_training.md)
  show map properties, component planes and the training progress.
- [`eemsom_cluster()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md)
  and
  [`eemsom_kselect()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_kselect.md)
  group the units by hierarchical clustering or k-means;
  [`eemsom_plot_silhouette()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_silhouette.md)
  shows silhouette widths.
- [`eemsom_split()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_split.md)
  splits an eemlist into groups for separate PARAFAC models;
  [`eemsom_predict()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_predict.md)
  maps new samples.
- [`eemsom_regions()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)
  groups wavelength pairs with correlated component planes into
  fluorescence regions (Ejarque-Gonzalez and Butturini 2014),
  [`ggeem()`](https://rdrr.io/pkg/staRdom/man/ggeem.html) plots them and
  [`eemsom_fri()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_fri.md)
  calculates area-normalised regional volumes (Chen et al. 2003).
- [`eemsom_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)
  assesses the stability of quantisation and neighbourhood relations by
  bootstrap or leave-one-out (de Bodt et al. 2002),
  [`eemsom_plot_stability()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_stability.md)
  compares it to a random map.
  [`eemsom()`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)
  stores its settings for this purpose.
- [`eempf_som()`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md)
  calculates a self-organising map from the Fmax loadings of a PARAFAC
  model (relative or absolute). Prototype EEMs are reconstructed from
  the PARAFAC components;
  [`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md)
  shows component planes of PARAFAC components with the new argument
  `comps`.
- The squared Euclidean distance is always used. kohonen otherwise
  switches to the Tanimoto distance for data in \[0, 1\] with row sums
  of 1.
- Vignette “Self-organising maps of EEM data” with the complete workflow
  from staRdom data to the interpretation of the map.
- [`eemsom_plot_cplane()`](https://matthiaspucher.github.io/staRsom/reference/eemsom_plot_cplane.md)
  gains `scale = "each"` to scale each component plane from 0 to 1.
- `citation("staRsom")` gives the package and the method paper of
  Ejarque-Gonzalez and Butturini (2014).

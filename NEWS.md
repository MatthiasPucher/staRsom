# staRsom 0.1.0

First release.

* `eem2mat()` and `mat2eem()` convert eemlists to matrices and back, including normalisation and handling of missing values.
* `eemsom()` calculates a self-organising map from EEM data; `eemsom_samples()`, `eemsom_codes()` and `eemsom_quality()` extract results.
* `ggeem()` method for `eemsom` objects plots the prototype EEMs arranged like the map; `eemsom_plot_map()`, `eemsom_plot_cplane()` and `eemsom_plot_training()` show map properties, component planes and the training progress.
* `eemsom_cluster()` and `eemsom_kselect()` group the units by hierarchical clustering or k-means; `eemsom_plot_silhouette()` shows silhouette widths.
* `eemsom_split()` splits an eemlist into groups for separate PARAFAC models; `eemsom_predict()` maps new samples.
* `eemsom_regions()` groups wavelength pairs with correlated component planes into fluorescence regions (Ejarque-Gonzalez and Butturini 2014), `ggeem()` plots them and `eemsom_fri()` calculates area-normalised regional volumes (Chen et al. 2003).
* `eemsom_stability()` assesses the stability of quantisation and neighbourhood relations by bootstrap or leave-one-out (de Bodt et al. 2002), `eemsom_plot_stability()` compares it to a random map. `eemsom()` stores its settings for this purpose.
* `eempf_som()` calculates a self-organising map from the Fmax loadings of a PARAFAC model (relative or absolute). Prototype EEMs are reconstructed from the PARAFAC components; `eemsom_plot_cplane()` shows component planes of PARAFAC components with the new argument `comps`.
* The squared Euclidean distance is always used. kohonen otherwise switches to the Tanimoto distance for data in [0, 1] with row sums of 1.
* Vignette "Self-organising maps of EEM data" with the complete workflow from staRdom data to the interpretation of the map.
* `eemsom_plot_cplane()` gains `scale = "each"` to scale each component plane from 0 to 1.
* `citation("staRsom")` gives the package and the method paper of Ejarque-Gonzalez and Butturini (2014).

# staRsom

Self-organising maps (SOM) for fluorescence excitation emission matrices (EEMs) of dissolved organic matter (DOM). staRsom is an add-on to [staRdom](https://cran.r-project.org/package=staRdom) and works with EEMs prepared there (`eemlist`).

staRsom follows the approach of Ejarque-Gonzalez and Butturini (2014):

* EEMs are unfolded, normalised and arranged on a SOM (`eemsom()`), calculated by the package [kohonen](https://cran.r-project.org/package=kohonen).
* The prototypes of the units are clustered to find groups of samples (`eemsom_kselect()`, `eemsom_cluster()`).
* Wavelength pairs with correlated component planes are grouped into fluorescence regions (`eemsom_regions()`) and quantified per sample (`eemsom_fri()`).
* The stability of the map is assessed by bootstrap or leave-one-out (`eemsom_stability()`).
* Groups of samples can be analysed separately by PARAFAC (`eemsom_split()`), new samples can be mapped (`eemsom_predict()`).
* Alternatively, a SOM is calculated from the loadings of a PARAFAC model (`eempf_som()`).

Plots: `ggeem()` for prototype EEMs and regions, `eemsom_plot_map()`, `eemsom_plot_cplane()`, `eemsom_plot_silhouette()`, `eemsom_plot_stability()`, `eemsom_plot_training()`.

## Installation

staRsom is not on CRAN yet. Install it from the source directory with

```r
# install.packages("remotes")
remotes::install_local("path/to/staRsom")
```

## Example

```r
library(staRsom)
data(eem_list)

eem_list <- eem_list |>
  eem_rem_scat(remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15, 15, 18, 19)) |>
  eem_interp(type = 1, extend = FALSE)

som <- eemsom(eem_list, seed = 42)
eemsom_plot_map(som, type = "umatrix")
ggeem(som)

som <- eemsom_cluster(som, k = 3)
eem_groups <- eemsom_split(eem_list, som, min_samples = 5)
```

The vignette shows the complete workflow:

```r
vignette("SOM_analysis_of_EEM", package = "staRsom")
```

## Citation

Please cite staRsom and the method paper of Ejarque-Gonzalez and Butturini (2014), see `citation("staRsom")`.

## References

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. https://doi.org/10.1371/journal.pone.0099618

Pucher, M., Wünsch, U., Weigelhofer, G., Murphy, K., Hein, T., Graeber, D. (2019). staRdom: Versatile Software for Analyzing Spectroscopic Data of Dissolved Organic Matter in R. Water 11(11), 2366. https://doi.org/10.3390/w11112366

Wehrens, R., Kruisselbrink, J. (2018). Flexible Self-Organizing Maps in kohonen 3.0. Journal of Statistical Software 87(7), 1-18. https://doi.org/10.18637/jss.v087.i07

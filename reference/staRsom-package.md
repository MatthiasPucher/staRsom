# staRsom: Self-Organising Maps for EEMs from DOM

Self-organising maps (SOM, Kohonen 2013)
[doi:10.1016/j.neunet.2012.09.018](https://doi.org/10.1016/j.neunet.2012.09.018)
for excitation emission matrix (EEM) data of dissolved organic matter
(DOM) prepared with 'staRdom'. Samples are arranged on a map, grouped by
clustering the prototype EEMs and fluorescence regions are derived from
correlated component planes following Ejarque-Gonzalez and Butturini
(2014)
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)
. Results can be split into sample groups for subsequent PARAFAC
analyses.

## Details

staRsom follows the approach of Ejarque-Gonzalez and Butturini (2014):
EEMs are arranged on a self-organising map
([`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md)),
the prototypes are clustered to find groups of samples
([`eemsom_cluster`](https://matthiaspucher.github.io/staRsom/reference/eemsom_cluster.md))
and wavelength pairs with correlated component planes are grouped into
fluorescence regions
([`eemsom_regions`](https://matthiaspucher.github.io/staRsom/reference/eemsom_regions.md)).
Please see the vignette for the complete workflow:
[`vignette("SOM_analysis_of_EEM", package = "staRsom")`](https://matthiaspucher.github.io/staRsom/articles/SOM_analysis_of_EEM.md).

## References

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

Pucher, M., Wünsch, U., Weigelhofer, G., Murphy, K., Hein, T., Graeber,
D. (2019). staRdom: Versatile Software for Analyzing Spectroscopic Data
of Dissolved Organic Matter in R. Water 11(11), 2366.
[doi:10.3390/w11112366](https://doi.org/10.3390/w11112366)

## See also

Useful links:

- <https://matthiaspucher.github.io/staRsom/>

- <https://github.com/MatthiasPucher/staRsom>

- Report bugs at <https://github.com/MatthiasPucher/staRsom/issues>

## Author

**Maintainer**: Matthias Pucher <matthias.pucher@ecosight.at>
([ORCID](https://orcid.org/0000-0002-9093-2199))

Authors:

- Matthias Pucher <matthias.pucher@ecosight.at>
  ([ORCID](https://orcid.org/0000-0002-9093-2199))

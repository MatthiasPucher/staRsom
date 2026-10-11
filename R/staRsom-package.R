#' @details staRsom follows the approach of Ejarque-Gonzalez and Butturini (2014): EEMs are arranged on a self-organising map (\code{\link{eemsom}}), the prototypes are clustered to find groups of samples (\code{\link{eemsom_cluster}}) and wavelength pairs with correlated component planes are grouped into fluorescence regions (\code{\link{eemsom_regions}}). Please see the vignette for the complete workflow: \code{vignette("SOM_analysis_of_EEM", package = "staRsom")}.
#'
#' @references Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and Correlation Analysis as a Tool to Explore Patterns in Excitation-Emission Matrix Data Sets and to Discriminate Dissolved Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618. \doi{10.1371/journal.pone.0099618}
#'
#' Pucher, M., Wünsch, U., Weigelhofer, G., Murphy, K., Hein, T., Graeber, D. (2019). staRdom: Versatile Software for Analyzing Spectroscopic Data of Dissolved Organic Matter in R. Water 11(11), 2366. \doi{10.3390/w11112366}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom eemR eem_names
#' @importFrom staRdom ggeem
## usethis namespace: end
NULL

# Calculate a self-organising map from EEM data

The EEMs are unfolded into vectors (see
[`eem2mat`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md)),
normalised and arranged on a self-organising map (SOM, Kohonen 2013).
Each unit of the map gets a prototype EEM (codebook vector). Similar
samples are mapped to the same or neighbouring units. The calculation is
done with [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html), please
see details there.

## Usage

``` r
eemsom(
  eem_list,
  dim = "auto",
  topo = "hexagonal",
  normalise = "sd",
  na = "remove",
  rlen = 500,
  alpha = c(0.05, 0.01),
  mode = "online",
  init = "pca",
  seed = NULL,
  cores = 1,
  verbose = FALSE,
  ...
)

# S3 method for class 'eemsom'
print(x, ...)
```

## Arguments

- eem_list:

  object of class eemlist. All EEMs need to have the same excitation and
  emission wavelengths and unique sample names.

- dim:

  size of the map, either `"auto"` (default) or a numeric vector
  `c(xdim, ydim)` with the number of columns and rows of units

- topo:

  topology of the map, `"hexagonal"` (default) or `"rectangular"`

- normalise:

  normalisation of each sample, see
  [`eem2mat`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md)

- na:

  handling of missing values, see
  [`eem2mat`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md)

- rlen:

  number of times the complete data set is presented to the map during
  training

- alpha:

  learning rate, vector of start and end value

- mode:

  training algorithm: `"online"` (default), `"batch"` or `"pbatch"`
  (parallel batch)

- init:

  initialisation of the codebook vectors: `"pca"` (default) places them
  on the plane of the first two principal components, `"random"` uses
  randomly selected samples

- seed:

  seed for the random number generator. With `NULL` (default), the
  current random number stream of the session is used, so
  [`set.seed`](https://rdrr.io/r/base/Random.html) before the call gives
  reproducible results. With a number, the seed is set for this
  calculation only and the random number state of the session is
  restored afterwards.

- cores:

  number of cores used with `mode = "pbatch"`

- verbose:

  print infos

- ...:

  additional parameters that are passed on to
  [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html), e.g. `radius`.
  `keep.data` is always TRUE, because the data are needed by
  [`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md)
  and other functions.

- x:

  object of class eemsom

## Value

object of class eemsom, a list containing

- som: object of class kohonen, see
  [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html)

- ex, em: excitation and emission wavelengths

- keep: logical vector of the unfolded wavelength pairs used in the SOM

- samples: sample names

- normalise, norm_factors: normalisation method and factors

- na: handling of missing values

- clusters: cluster of each unit, NULL until the units are clustered

- settings: arguments used for
  [`som`](https://rdrr.io/pkg/kohonen/man/supersom.html), needed by
  [`eemsom_stability`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md)

## Details

With `dim = "auto"`, the number of units is \\5 \sqrt{n}\\ (Vesanto and
Alhoniemi 2000) but not more than the number of samples n. The ratio of
the side lengths follows the two largest eigenvalues of the data
(Bieroza et al. 2009, Ejarque-Gonzalez and Butturini 2014): like in the
SOM Toolbox (Vesanto et al. 2000), it is the square root of the ratio of
the eigenvalues, i.e. the ratio of the standard deviations along the
first two principal components. Ratios above the square root of the
number of units are replaced by 1 to avoid very elongated maps.

The default learning rate (0.05 to 0.01) and the default neighbourhood
radius (decreasing from two thirds of all unit distances to zero) of
[`som`](https://rdrr.io/pkg/kohonen/man/supersom.html) correspond to the
settings of Ejarque-Gonzalez and Butturini (2014).

The squared Euclidean distance (`dist.fcts = "sumofsquares"`) is always
used. Otherwise, kohonen would switch to the Tanimoto distance for data
between 0 and 1 with row sums of 1, which occurs with
`normalise = "sum"`. With `na = "keep"`, distances are calculated from
the available wavelength pairs of a sample and scaled up to the number
of all wavelength pairs.

The initialisation with `init = "pca"` is deterministic. The order of
samples in the online training is random, so please set `seed` or use
[`set.seed`](https://rdrr.io/r/base/Random.html) before the call for
reproducible results.

The training time grows with the numbers of samples, units and
wavelength pairs. As a rough guide for one core, online training with
`rlen = 500` takes about a minute for 200 EEMs with 3000 wavelength
pairs and several minutes for 1000 EEMs. `mode = "pbatch"` with
`cores > 1`, a smaller `rlen` or a lower resolution of the EEMs (see
[`eem_range`](https://rdrr.io/pkg/staRdom/man/eem_range.html) and
[`eem_exclude`](https://rdrr.io/pkg/staRdom/man/eem_exclude.html);
Ejarque-Gonzalez and Butturini (2014) used 366 wavelength pairs) reduce
the time. This matters for
[`eemsom_stability`](https://matthiaspucher.github.io/staRsom/reference/eemsom_stability.md),
which repeats the training many times.

## References

Kohonen, T. (2013). Essentials of the self-organizing map. Neural
Networks 37, 52-65.
[doi:10.1016/j.neunet.2012.09.018](https://doi.org/10.1016/j.neunet.2012.09.018)

Vesanto, J., Alhoniemi, E. (2000). Clustering of the self-organizing
map. IEEE Transactions on Neural Networks 11(3), 586-600.
[doi:10.1109/72.846731](https://doi.org/10.1109/72.846731)

Vesanto, J., Himberg, J., Alhoniemi, E., Parhankangas, J. (2000). SOM
Toolbox for Matlab 5. Report A57, Helsinki University of Technology.

Bieroza, M., Baker, A., Bridgeman, J. (2009). Exploratory analysis of
excitation-emission matrix fluorescence spectra with self-organizing
maps as a basis for determination of organic matter removal efficiency
at water treatment works. Journal of Geophysical Research:
Biogeosciences 114(G4), 2009JG000940.
[doi:10.1029/2009JG000940](https://doi.org/10.1029/2009JG000940)

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

Wehrens, R., Kruisselbrink, J. (2018). Flexible Self-Organizing Maps in
kohonen 3.0. Journal of Statistical Software 87(7), 1-18.
[doi:10.18637/jss.v087.i07](https://doi.org/10.18637/jss.v087.i07)

## See also

[`som`](https://rdrr.io/pkg/kohonen/man/supersom.html),
[`eemsom_samples`](https://matthiaspucher.github.io/staRsom/reference/eemsom_samples.md),
[`eemsom_codes`](https://matthiaspucher.github.io/staRsom/reference/eemsom_codes.md),
[`eemsom_quality`](https://matthiaspucher.github.io/staRsom/reference/eemsom_quality.md)

Other self-organising maps:
[`eem2mat()`](https://matthiaspucher.github.io/staRsom/reference/eem2mat.md),
[`eempf_som()`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md),
[`mat2eem()`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md)

## Examples

``` r
data(eem_list)

som <- eemsom(eem_list, rlen = 100, seed = 42)
som
#> Self-organising map of 15 EEMs
#>  map: 3 x 5 hexagonal units
#>  wavelength pairs: 4554 of 4554 
#>  normalisation: sd 
#>  clusters: none 
eemsom_samples(som)
#>          sample unit   x         y  distance hits  umatrix cluster
#> 1      d0680sfK    9 3.5 2.5980762  2.300073    1 36.48733      NA
#> 2       d1266sf    3 3.5 0.8660254  9.054573    2 10.92376      NA
#> 3      d1268sfK    3 3.5 0.8660254 11.834547    2 10.92376      NA
#> 4      d1543sfK   13 1.5 4.3301270  7.938452    5 11.58416      NA
#> 5      d1549sfK   13 1.5 4.3301270  9.879895    5 11.58416      NA
#> 6      d1550sfK    7 1.5 2.5980762  4.606082    3 15.87562      NA
#> 7      d1557sfK   13 1.5 4.3301270  4.190952    5 11.58416      NA
#> 8      d1559sfK   13 1.5 4.3301270  8.586989    5 11.58416      NA
#> 9      d1565sfK    7 1.5 2.5980762  6.557108    3 15.87562      NA
#> 10     d1567sfK   13 1.5 4.3301270 18.078294    5 11.58416      NA
#> 11     d1568sfK    7 1.5 2.5980762  8.106521    3 15.87562      NA
#> 12       d492sf    1 1.5 0.8660254  3.276772    2 21.19327      NA
#> 13       d661sf    1 1.5 0.8660254  3.092985    2 21.19327      NA
#> 14       d667sf    2 2.5 0.8660254  1.694963    1 14.83993      NA
#> 15 dSFB445pwd1t   15 3.5 4.3301270  2.842423    1 55.31539      NA
```

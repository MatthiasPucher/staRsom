# Unfold EEMs into a matrix of samples and wavelength pairs

Each EEM in an eemlist is unfolded into a vector. The result is a matrix
with one row per sample and one column per pair of excitation and
emission wavelengths. This matrix is the input for the self-organising
map in
[`eemsom`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md).
The unfolding can be reversed by
[`mat2eem`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md).

## Usage

``` r
eem2mat(eem_list, normalise = "sd", na = "remove")
```

## Arguments

- eem_list:

  object of class eemlist. All EEMs need to have the same excitation and
  emission wavelengths and unique sample names.

- normalise:

  normalisation of each sample: `"sd"` (default) divides each EEM by its
  standard deviation like
  [`norm_array`](https://rdrr.io/pkg/staRdom/man/norm_array.html),
  `"max"` by its maximum (Ejarque-Gonzalez and Butturini 2014), `"sum"`
  by the sum of its fluorescence and `"none"` keeps the original values.

- na:

  handling of missing values (e.g. removed scatter): `"remove"`
  (default) removes wavelength pairs that are missing in at least one
  sample, `"keep"` keeps wavelength pairs that are missing in some of
  the samples. Wavelength pairs that are missing in all samples are
  always removed, see details.

## Value

A matrix with samples in rows and wavelength pairs in columns.
Attributes: `ex` and `em` (wavelengths of the original EEMs), `keep`
(logical vector stating which of the unfolded wavelength pairs are
columns of the matrix), `samples`, `normalise` and `norm_factors`.

## Details

Columns are named "ex\|em", e.g. "250\|450". The EEM matrix is unfolded
column-wise, so emission wavelengths change faster than excitation
wavelengths.

Normalisation removes differences in concentration and focuses on the
shape of the EEMs. The normalisation factors are calculated from the
remaining wavelength pairs (after removing missing values) and are saved
as an attribute. They need to be positive, so EEMs with a maximum or a
sum of zero or below cannot be normalised with `"max"` or `"sum"`.

Wavelength pairs that are missing in all samples (e.g. scatter that was
removed but not interpolated) carry no information and are removed even
with `na = "keep"`. Otherwise, the self-organising map would contain
arbitrary values for them, because kohonen fills missing values in the
initial codebook vectors with random values. With `na = "keep"`, kohonen
calculates the distances between samples and codebook vectors from the
available wavelength pairs and scales them up to the number of all
wavelength pairs.

## References

Ejarque-Gonzalez, E., Butturini, A. (2014). Self-Organising Maps and
Correlation Analysis as a Tool to Explore Patterns in
Excitation-Emission Matrix Data Sets and to Discriminate Dissolved
Organic Matter Fluorescence Components. PLoS ONE 9(6), e99618.
[doi:10.1371/journal.pone.0099618](https://doi.org/10.1371/journal.pone.0099618)

## See also

[`mat2eem`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md),
[`eem2array`](https://rdrr.io/pkg/staRdom/man/eem2array.html)

Other self-organising maps:
[`eempf_som()`](https://matthiaspucher.github.io/staRsom/reference/eempf_som.md),
[`eemsom()`](https://matthiaspucher.github.io/staRsom/reference/eemsom.md),
[`mat2eem()`](https://matthiaspucher.github.io/staRsom/reference/mat2eem.md)

## Examples

``` r
data(eem_list)

m <- eem2mat(eem_list)
dim(m)
#> [1]   15 4554
attr(m, "norm_factors")
#>     d0680sfK      d1266sf     d1268sfK     d1543sfK     d1549sfK     d1550sfK 
#>   0.08683636   0.08967149   0.09703871   0.06498094   0.04673046   0.05382821 
#>     d1557sfK     d1559sfK     d1565sfK     d1567sfK     d1568sfK       d492sf 
#>   0.05307236   0.03873703   0.04166047   0.09772790   0.05926670   0.06871059 
#>       d661sf       d667sf dSFB445pwd1t 
#>   0.06043481   0.07290243   0.11637892 
```

---
title: NEWS
editor_options:
  markdown:
    wrap: 72
---

# photobiology 0.13.2

* Support `log()` and `sqrt()` transformations for `local.reference` in
`find_peaks()` and `find_valleys()`.
* Add `s_quantile()` methods for parallel quantiles from a collection of
spectra.
* Fix long-standing bug in interpolation of spectra: `approx()` called instead
of `spline()` even for sparse or short input. (Values returned by 
`interpolate_wl()` and low-level functions it calls are different in
some cases.)

# photobiology 0.13.1

* Rewrite function `find_peaks()` fixing a bug in the logic of threshold tests.
This modifies the behaviour compared to version 0.13.0, when first introduced.
* Add function `find_valleys()` and edit `valleys()` methods to use it.
* Add utility function `check_wl_stepsize()`.
* These changes, repair the behaviour of `peaks()` and `valleys()` methods. In
particular, they remain backwards compatible with versions < 0.13.0.
* Update the vignette.

# photobiology 0.13.0

This update focus is on 1) more efficient storage of metadata in attributes, 2)
an improved match between operations on spectra stored as collections and
multiple spectra stored in long form, and 3) finer control of the detection and
extraction of peaks and valleys. It also includes a bug fix resolving a
problem in package 'ooacquire'. R (>= 4.1.0) is now required.

* Support simplification of attributes in `what_measured()`, `when_measured()`, 
`where_measured()`, and `how_measured()` attribute accessors.
* Improve copying of attributes in row-wise summary methods: `"when.measured"`,
 `"where.measured"`, `"what.measured"`, and `"how.measured"` attributes are
copied, and when unique across the summarised spectra, they are simplified. This
change affects specializations for collections of spectra of methods `s_mean()`,
`s_median()`, `s_prod()`, `s_range()`, `s_sd()`, `s_se()`, `s_sum()`. `s_var()`,
`s_mean_se()`.
* Simplify repeated identical attributes when printing multiple spectra in long 
form or as a collection, and when printing their summaries. Simplification
controlled by formal parameter `attr.simplify` with default to `TRUE`.
* Handle correctly simplified attributes in query and print methods.
* Accept spectra in long form in addition to collections of spectra as input to
`add_attr2tb()` and `spct_metadata()`. 
* Methods `peaks()`, `valleys()`, `wls_at_target()`, `spikes()`, and `despike()`
warn when called with spectral data expressed on irregular wavelength steps,
such as after applying `thin_wl()`.
* Implement a local (within-window span) threshold for peak height and valley 
depth in `stat_peaks()`  and `stat_valleys()`, using parameters 
`local.threshold` and `local.reference`.
* **Code breaking:** Rename parameter `ignore_threshold` into `global.threshold`
in `find_peaks()`, `get_peaks()`, `peaks()`, `get_valleys()`  and `valleys()`
for naming consistency and clarity.
* The scaling applied to user-supplied values for `global.threshold` and 
`local.threshold` can be controlled by passing a `character` argument to 
`threshold.scaling`. Non-scaled thresholds are also supported.
* Bug fix: in `trimInstrDesc()` and `trimInstrSettings()` handle correctly
missing record fields.

# photobiology 0.12.0

**This is a major update with new features and includes significant changes to
internal code and the migration of some functions to other packages. As of
2024-12-26 the development version of 'photobiology' (>= 0.11.4.9004) is
compatible with the current CRAN versions of 'ggspectra' (>= 0.3.14) and of
other packages in the suite. Some functions have been migrated to package
'photobiologySun' and others to 'SunCalcMeeus'.**

- Redesign the user interface for normalization and conversions of units or 
quantities. These changes keep the logic of both package code and user code 
simpler. These changes can break some user code.
  * `normalize()` methods no longer supports on-the-fly change of units or 
  quantities.
  * Methods for conversion between quantities and between units update the
  existing normalization.
  * Support in `normalize()` the undoing of an existing normalization with 
`norm = "undo"`. Available when normalization has been done with 'photobiology'
(>= 0.10.9), as earlier versions did not store the normalization multipliers as 
metadata.
  * Support normalization of multiple columns in spectra, such as photon- and
energy irradiances.
  * Query method `normalization()` implemented for spectra, summaries of spectra
and collections of spectra as equivalent to `getNormalization()`.
- Update range-checks of spectra to tolerate 1 in 250 pixels "mildly" off-range,
except `cps_spct` where 1 in 100 is tolerated and `raw_spct` with no range check
for counts. Update the messages to report the number of off-range values in
addition to the extreme values.
- Add support for `attenuation.mode = "scattering"` in filter properties 
attribute.
- Move to package ['SunCalcMeeus'](https://docs.r4photobiology.info/SunCalcMeeus/) 
the functions and methods for sun position and day length calculations.
- Move functions related to energy, water and carbon exchange between
vegetation and the atmosphere, including those for evapotranspiration rates,
water content in the atmosphere and radiation balance functions to package
'photobiologySun'. This is a code breaking change that will require adding
`library(photobiologySun)` to scripts that call these functions.

# photobiology 0.11.4

- Fix bug in function `subset2mspct()`: failure with same value of attribute
`when.measured` in all spectra.
- Fix bug in function `rbindspct()`: unnecessary setting of missing metadata
to `NA`s.
- Fix bug in `normalize()` methods: failure of repeated updates.
- Implement methods `where_measured()`, `what_measured()`, `how_measured()` and
`when_measured()` for data frames (and tibbles).
- Implement previously missing `setHowMeasured()` and `setWhatMeasured()` 
methods for collections of spectra.
- Support negative wavelengths as input in `color_of.numeric()` by changing 
sign. This is temporary patch to allow reverse transform in scale in 'ggspectra'.
- Replace use of unsupported time zone name "EET" by "Europe/Helsinki" for
compatibility with future R (==4.5.0).
- Update CIE data to the most recent version at 1 nm for A and D65 illuminants
and add the D50 illuminant.
- Rebuild all data objects.

# photobiology 0.11.3

Bug fixes and improved printing of spectra and their summaries. This changes
the printed header text, but should not break code.

- Fix bug that prevented use of `s_mean_se_band()`.
- Implement methods `s_mean()`, `s_mean_se()`, `s_mean_se_band()`, 
`s_median()`, `s_sd()`, `s_se()`, `s_var()`, `s_sum()`, `s_prod()`, `s_range()` 
for class `generic_spct`.
- Support renaming of the `idfactor` with `setIdFactor()`.
- New method `make_var_labels()` dynamically creates a named list of labels for 
the variables in spectral objects [labels' texts may change in next version].
- Update `print.generic_spct()` to include variable labels in header, instead
of showing metadata directly.
- Update `summary.generic_spct()` to save variable labels in returned object.
- Update `print.summary_generic_spct()` to print a header equal to the one
printed by the updated `print()` methods.
-  Update `summary.generic_spct()` to optionally expand objects containing
multiple spectra in long form into collections of spectra in advance of creating
the summary.
-  Update `summary.generic_mspct()` to optionally expand member objects and
call `summary.generic_spct()` on each of them.

# photobiology 0.11.2

Mostly a bug-fix release.

- Implement missing methods `normalize.solute_mspct()` and 
`normalize.generic_mspct()` 
- Accept arbitrary function names in addition to function objects in 
`fscale()` methods.
- Fix major bug affecting `fscale()` methods when applied to `<xxx>_spct` objects
containing multiple spectra in long form.
- Fix bug in `rbindspct()` non-default `idfactor` values could be ignored in
some cases.
- Fix bug in `normalise()`, `fscale()`, and `smooth_spct()` methods: dropped 
`idfactor` attribute when operating on `<xxx>_spct` objects containing multiple 
spectra in long form.
- Fix bug triggered in R < 4.3.0: missing argument for `origin` in two calls to
`as.POSIXct()`.

# photobiology 0.11.1

The main enhancements in this update are **1)** the implementation of proper
handling of metadata attributes in objects containing multiple spectra in long
form and during their two way conversion to collections of spectra, and **2)**
improved performance of the computation of irradiances for spectra stored in
long or tidy form.

- Implement subsetting of metadata in `subset()` method and the extraction 
operator `[]` for objects of class `generic_spct` and derived classes. This 
makes it possible to extract a subset of spectra from an object containing
multiple spectra in long form, such as a time series of spectra, with no
extraneous metadata attribute values carried along. In earlier 
versions metadata were not subset.
- In `rbindspct()` implement simplification of metadata by removing unnecessary
duplication of invariant values.
- Rewrite method `pull_sample()` specialization for `generic_spct` for faster 
performance.
- Revise `irrad()`, `e_irrad()` and `q_irrad()` adding parameter `return.tb`
making it possible to force the return of a tibble even for summaries of
individual spectra.
- Revise `irrad()`, `e_irrad()` and `q_irrad()` for faster performance with
multiple spectra in long-form. (Adding attributes to the returned tibble is not 
yet supported except for `when.measured`.)
- Revise `normalise()`, `fscale()` and `smooth_spct()` methods to accept
multiple spectra in long form as their argument returning a similar object, with
the spectra individually normalised or scaled.
- Revise `peaks()`, `valleys()`, `wls_at_target()`, `spikes()` and `despike()`
methods to consistently return an object as the same class as their `x` 
argument. 
- Revise `when_measured2tb()` to support `generic_spct` in addition to 
`generic_mspct` objects. (Methods to add other attributes are not yet revised.)
- Add parameter `span` to `thin_wl()` methods, with the previously hard-coded
value of 21 as default.
- Improve handling of `idFactor` and other non-numeric variables in `clip_wl()`
and `trim_wl()` so that they are filled with good values instead of NA when 
possible and fix a bug that converted non-numeric variables into numeric ones
in the returned value.
- Revise `getNormalized()`, `getNormalised()`, `getNormalization()`, 
`getNormalisation()`, `is_normalized()`, and `is_normalised()` to support
collections of spectra as their argument in addition to individual spectra. In
this case they return a named list.
- Add a short time series of sunlight spectra in objects `sun_evening.spct` and
`sun_evening.mspct`.
- When checking `raw_spct` and `cps_spct` do not emit a message about renaming 
columns if option `photobiology.verbose` is set to `FALSE`.
- Fix bug: some operations failed to copy all metadata attributes to the
returned value. Some of the dropped attributes are used in package 'ooacquire' 
(>= 0.4.1).
- Fix bug affecting `irrad()`, `e_irrad()` and `q_irrad()` causing a crash with
argument `use.hinges = FALSE`.
- Fix bug affecting handling of existing normalizations in conversions between
quantities used for describing spectral properties filters, visible only as a
spurious warning when plotting with 'ggspectra'.
- Fix handling of default `idFactor` in `subset2mspct()`.
- Spectral data objects included in the package and used in examples and unit
tests have been rebuilt making small changes in textual metadata, and/or adding
new fields with additional information.
- Values and length of objects `polyester.spct` and `yellow_gel.spct` have
changed slightly as storage use was optimized. This is a code breaking change
if positional indexes have been used, but otherwise differences in computed
values are minor.

# photobiology 0.11.0

- Improve handling of objects containing multiple spectra in long form. Most 
methods automatically expand these objects into collections of spectra and
recursively dispatch themselves. Many operations that would earlier fail with
an error message are now handled transparently.
- Include 'entrance.optics' in the output of `print.instr_desc()`.
- Update`trimInstrDesc()` so that by default it keeps the `entrance.optics`
  field.
- Add method `pull_sample()` for pulling random samples of spectra, replacing
functions functions `sample_spct()` and `sample_mspct()` added in 0.10,17.
- In `normalized_diff_ind()` methods, rename `plus.w.band` into `w.band.plus` 
and `minus.w.band` into `w.band.minus` for consistency with `Tfr_normdiff()`,
`Rfr_normdiff()` and other functions.
- Spectral data objects included in the package and used in examples and unit
tests have been rebuilt leading to small changes in data values. New data have
been added to demonstrate and test recently added features.

# photobiology 0.10.17

- Fix bug in `subset2mspct()` when applied to the case of a collection
  containing a single spct object containing a single spectrum.
- Add functions `sample_spct()` and `sample_mspct()`.
- Add specialization of `smooth_spct()` method for class `"cps_spct"`.
- Add `e_fraction()` and `q_fraction` methods.
- Add `Rfr_ratio()`, `Rfr_fraction()` and `Rfr_normdiff()`.
- Add `Tfr_ratio()`, `Tfr_fraction()` and `Tfr_normdiff()`.
- Add formal parameter `quantity` to `q_ratio()` and `e_ratio()` methods.
- Extend the "update" of normalizations to objects of class 
  `"filter_spct"` and fix a bug affecting setting of metadata attribute during
 updates to the normalization of `"source_spct"` and `"response_spct"`
  objects.

# photobiology 0.10.16

- Fix bug in `subset2mspct()` introduced in 0.10.15 affecting collections of 
  spectra with a single spct object with multiple spectra in long form as only 
  member.
- Fix failure of automatic registration of methods `log2()` and `log10()`.

# photobiology 0.10.15

- Add helper function `spct_wide2long()` a simple pure R replacement for
  `dplyr::pivot_longer()` for spectra.
- Revise function `subset2mspct()` to accept collections of spectra and
  subset (split) members containing multiple spectra in long form.
- Revise function `subset2mspct()` to accept any valid `spct` object and
  always return a `mspct` object.
- Add function `is_daytime()` a wrapper on `day_night()` returning a logical
  vector.
- Make the name and label returned by a call to `waveband()` accepting all
  defaults argument more informative.
- Update `e2q()` and `q2e()` to better handle previously normalized spectra 
  passed as arguments: re-normalization is applied by default.
- Increase accuracy of conversions in `e2q()` and `q2e()`. The change in 
  computed values is at most 20 parts per million.
- Add functions for conversions among quantities representable as wave length:
  `wl2wavenumber()`, `wavenumber2wl()`, `wl2frequency()`, `frequency2wl()`,
  `wl2energy()` and `energy2wl()`.

# photobiology 0.10.14

- Fix several bugs created by code-breaking changes in 'tidyselect' 1.2.0, and
  possibly 'rlang' 1.0.6, affecting 'dplyr'.
- Fix other bugs, including handling of spectra with no non-missing data.
- Add function `illuminance()`.

# photobiology 0.10.13

- Improve handling of missing and default wavebands in `irrad()` and in
  `trim_waveband()`.
- Add preliminary support for filter stacks in `filter.properties` attribute.
- Add `summary()` method for collections of spectra.
- **Fix bug** in extract (`[`) operator for collections of spectra, resulting in
wrong values for dimension attribute (`"mspct.dim"`).

# photobiology 0.10.12

------------------------------------------------------------------------

Conversions between `Date` and `POSIXct` objects are tricky because the
former do not store information on the time zone. A change in
'lubridate' 1.8.0 made a previously working approach to these
conversions silently fail to apply the shift to the hours. *In the
current version of 'photobiology', if no time zone argument is passed
concurrently with a date, the date is assumed to be in UTC. If this time
zone does not match the location given by the geocode, the date used for
the calculations can be wrong by one day.*

------------------------------------------------------------------------

This release corrects problems triggered by recent updates to packages
'lubridate' and possibly 'tibble' (reported by *putmanlab* in [issue
#7](https://github.com/aphalo/photobiology/issues/7)) and adds
enhancements for class `solute_spct`.

-   Bug fixed: With 'lubridate' (1.8.0) but not with previous versions,
    functions `day_night()`, `sunrise_time()`, `noon_time()` and
    `sunset_time()` would return wrong time values when non-default
    arguments to parameter `tz` were passed together with objects of
    class `Date` passed as arguments to `date`.
-   Add methods `as.filter_spct()` and `as.solute_spct` specialised for
    two-way conversion between objects of classes `solute_spct` and
    `filter_spct`.
-   Revise the class `solute.properties` adding fields `solvent.name`
    and `solvent.ID`.
-   Revise documentation checking that units expected for arguments and
    of returned values are clearly indicated and correctly formatted.
    Update outdated text and correct mistakes and revise unclear
    explanations.

# photobiology 0.10.11

-   Add new classes of objects `solute_spct` and `solute_mspct` to be
    used to store molar (default) and mass based coefficients of
    attenuation describing overall attenuation, or attenuation by
    absorption or by scattering. Implement the corresponding methods.
    (Unstable: interface may change).
-   Add example data for two substances: `water.spct` and
    `phenylalanine.spct`.
-   Rewrite `join_mspct()` to use interpolation when wavelengths differ
    among member spectra. This should not break old code but output can
    slightly differ.
-   Expand syntax accepted for `character` arguments passed to parameter
    `target` in all `wls_at_target()` methods.
-   Fix failure to handle spectra with zero rows, a bug affecting
    several methods, operators and functions including `rbindspct()` and
    `find_wls()`.
-   Fix bug in `rowwise_filter()` affecting parallel summaries of
    absorptance.
-   Fix bugs in extraction and replacement functions for collections of
    spectra, possibly triggered by changes in R \>= 4.0.0.
-   Add method `s_mean_se_band_band()`.

# photobiology 0.10.10

-   Update `normalize()` methods to support updating an already present
    normalization (`norm = "update"`) and skipping the normalization
    altogether (`norm = "skip"`).
-   Update `normalize()` methods to store `range` in the attribute, and
    `getNormalized()` to return it.
-   Update `normalize()` methods to correctly handle normalization of
    previously normalized spectra, and add flexibility to the
    normalization of previously scaled spectra.
-   Add `getScaling()` and fix minor inconsistency in value returned by
    `getScaled()`.
-   Fix bug in `getNormalization()` (wrong named member in returned
    value from spectra with no normalization data).
-   Fix bug resulting in `"normalization"` attribute not being copied.
-   Fix bug resulting in not all relevant attributes being copied to the
    value returned by `summary.generic_spct()`.
-   Improve printing of metadata for normalization and rescaling.
-   Fix bug in `shared_member_class()` (wrong value returned for empty
    collections).
-   Update `smooth_spct()` to handle bad arguments passed to `method`
    without crashing and add support for skipping smoothing
    (`method = "skip"`).

# photobiology 0.10.9

-   Update `smooth_spct()` methods so that `NA` values in `wl.range` are
    handled as documented and consistently with other methods in the
    package.
-   Update to accommodate code-breaking change in 'dplyr' (\>= 1.0.8).

# photobiology 0.10.8

-   Update functions `normalize()`, `setNormalized()` and
    `getNormalized()`, and add new function `getNormalization()`. These
    changes implement the storage in attribute `normalization` of the
    operation done.
-   Fix bug in `mat2mspct()` affecting matrices with more than 26
    columns and without `colnames` previously set.
-   Fix bug in `rowwise` methods.

# photobiology 0.10.7

-   Add function `ET_ref()` for computation of reference
    evapotranspiration, implementing the original FAO56 formulation of
    the Penman-Monteith method as well as modified in 2005 for tall and
    short vegetation according to ASCE-EWRI. The formulation is that for
    ET expressed in mm/h, but modified to use as input flux rates in
    W/m2 and pressures expressed in Pa.
-   Add function `net_radiation()` that computes the long wave net
    radiation balance if down-welling long wave radiation is available
    and otherwise estimates it.
-   Add function `irrad_extraterrestrial()` that computes down-welling
    solar irradiance on a horizontal plane at the top of the atmosphere.
-   Revise function `sun_angles()` to also return the Sun to Earth
    distance.

# photobiology 0.10.6

-   Fix boundary-case bug in `msmsply()`.
-   Revise the computation of the default for `dyn.range` in `cps2Tfr()`
    and `cps2Rfr()` so that it takes into account the relative signal in
    the reference spectrum.
-   Add parameter `missing.pixs` to `cps2irrad()` so that corrupted
    too-short spectra can be converted if the location of missing pixels
    is known.
-   Add row-wise summaries for `raw_mspct` and `cps_mspct` objects.
-   Add support of multiple spectra in long form to `irrad()`,
    `e_irrad()`, `q_irrad()`, `q_ratio()`, `e_ratio()`, `qe_ratio()`,
    `eq_ratio()`, `absorbance()`, `absorptance()`, `trasmittance()`,
    `reflectance()` methods.
-   Add warning for handling of multiple spectra in long form to
    `integrate_spct()` method.
-   Fix handling of `na.rm = TRUE` in `find_peaks()`.
-   **Move git repository from Bitbucket to Github.**
-   Set up Github action for CRAN-checks on Windows, OS X and Ubuntu.

# photobiology 0.10.5

-   Make the precomputed color data private.
-   Add method `drop_user_cols()` to remove user-defined columns from
    spectra.
-   Add method `collect2mspct()` and rename method `uncollect()` into
    `uncollect2spct()`.
-   Add method `spct_metadata()` to query the value of metadata
    attributes.
-   Revise `add_attr2tb()` expanding support to all metadata attributes.
-   Revise `rbindspct()` to gracefully handle duplicate member names in
    its input.
-   Revise `smooth_spct()` methods: new parameter wl.range and bug fix
    handling of strength in `"custom"` method.
-   Revise `compare_spct()` function to accept scaled and normalized
    spectra with a warning.
-   Implement attribute `"response.type"` to distinguish between
    response spectra and action spectra stored in `response_spct`
    objects. Add methods `setResponseType()` and `getResponseType()`.

# photobiology 0.10.4

-   Bug fixes.

-   Improved performance in color-related functions, mainly benefiting
    package 'ggspectra'.

-   Handle gracefully and consistently special input in
    `fast_color_of_wl()`.

-   Add `fast_wb2rect_spct()`, which uses precomputed color definitions
    for narrow wavebands and optionally simplifies the returned spectrum
    by merging neighboring rectangles of identical color.

-   Add `fast_color_of_wb()` that uses precomputed color definitions for
    narrow wavebands.

-   Add parameter force to `check_spct()` methods, so that critical
    checks cannot be disabled.

-   Implement math functions for class `generic_spct`.

-   BUG FIX: ERROR in CRAN check because of bad example in docs.

-   BUG FIX: `tag()` would fail to assign `wb.color` and `wb.name` to
    longest `w.length` value in spectrum.

# photobiology 0.10.3

-   Handle gracefully bad data input in `normalised_diff_ind()`.
-   Implement Fresnel's formulae for computation of reflectance of a
    plane interface from relative refractive index.
-   Implement Fraunhofer's formulae for computation of diffraction in a
    single slit and diffraction plus interference in a double slit.
-   Enhance `as_tod()` and implement `format()` and `print()` methods
    for *time-of-day*.
-   Update `tag()` to use precomputed color definitions, when possible,
    to improve performance.
-   BUG FIX: Remove bad class exports from NAMESPACE.
-   New features of dplyr (\>= 1.0.0) are used, so this new version is
    required.

# photobiology 0.10.2

-   Fix bug in `color_of()`.
-   Fix bug in `merge2object_spct()`.
-   Fix bug in `merge_attributes()`.
-   Fix bug in `interpolate_spct()`.
-   Fix bug in Extract `[ ]`.
-   Rebuild white LED example spectra.

# photobiology 0.10.1

-   This update brings improved handling of conversions among quantities
    used to describe filters. This required adding a mechanism to store
    properties of filters as metadata, revising the code behind
    mathematical operations and mathematical functions. Also adding
    conversion methods and obviously methods for setting and retrieving
    the new metadata. It also brings support for arbitrary chromaticity
    definitions to colour handling and tagging of spectra.
-   Reshuffling the code and testing brought to light several bugs and
    rough edges in infrequently used functions, which I have corrected.
-   Some methods gained new specializations for data frames, which were
    needed to streamline some of the statistics in package 'ggspectra'.
-   This update also brings some changes to error checking and messages,
    including changing several earlier warnings into errors.

### New

-   Add `setFilterProperties()`, `getFilterProperties()`,
    `filter_properties()` and `filter_properties<-()`.
-   Add `convertTfrType()`, `convertThickness()`.
-   Add `Afr2T()`, `any2T()`, `any2A()`, and `any2Afr()`.
-   Add `print()` method for filter properties.
-   Update example data for filters by adding filter properties.

### Fix bugs and "polish rough edges".

-   **Major bug** in `T2Afr()` was causing wrong values to be returned!!
-   **Major bug** in `clean.object_spct()`.
-   Rewrite much of the code for dispatching math operations with
    spectral objects as operands.
-   Revise `find_peaks()` so that `ignore.threshold` and `strict`
    arguments are respected also when `span = NULL`. This affects all
    peaks- and valleys-related methods. Changes return values for the
    previously undocumented use of negative threshold values and with
    `span = NULL` and `strict = TRUE` as arguments.
-   Add `wls_at_target()` method for data frames.
-   Revise `wls_at_target()` methods for consistency, add examples.
-   Revise `rbindspct()` to not add columns to returned value when input
    is consistent.
-   Revise `rbindspct()` to allow control of metadata copying.
-   Add `is_absorptance_based()`\`.

### Defunct

-   `T2T()`, `setAfrType()`, `getAfrType()`. No deprecation as not used.

------------------------------------------------------------------------

Full functionality of `convertThickness()` and `convertTfrType()`
requires that `filter_spct` objects have the new `"filter.properties"`
attribute or that they contain a column `"Rfr"` with spectral
reflectance data like `object_spct` objects do.

------------------------------------------------------------------------

# photobiology 0.10.0

### Enhance functionality of peak- and valley-related functions:

-   Revise `peaks()` and `valleys()` methods to support peak fitting.
    Currently only spline interpolation is implemented. Default
    behaviour remains unchanged. (The implementation of peak fitting is
    still experimental and returned values may not be reproducible with
    future versions of the package.)

### Enhancements to summary functions:

-   Revise `q_ratio()`, `e_ratio()`, `qe_ratio()` and `eq_ratio()`
    adding new parameter `scale.factor`.

-   Revise `irrad()`, `q_irrad()` and `e_irrad()`; `fluence()`,
    `q_fluence()`, and `e_fluence()`; `q_ratio()`, `e_ratio()`,
    `qe_ratio()` and `eq_ratio()`; `response()`, `q_response()`, and
    `e_response()` methods so that they add shorter but still
    informative names to returned numeric values and to columns in
    returned data frames. Add formal parameters naming and name.tag to
    allow user control of the names.

-   Revise `absorbance()`, `absorptance()`, `transmittance()` and
    `reflectance()` so that they create shorter, and more informative
    names for returned values.

-   Revise `msdply()` to treat as special cases all the methods
    described above, to support the especializations of these methods
    for collections of spectra. For those functions not handled as
    special cases, returned values are, as earlier, tagged with the name
    of the applied function.

    ### Add support for extraction of spikes and for despiking:

-   Add function `find_spikes()` for numeric vectors, useful for Raman
    spectra.

-   Add function `replace_bad_pixs()` operating on numeric vectors.

-   Add method `spikes()` for vectors, data frames and spectra.

-   Add method `despike()` for vectors, data frames and spectra.

    ### Add example data:

-   Add data for cone fundamentals for human vision.

    ### Fix bugs:

-   Revise `trim_waveband()` so that it preserves the names of list
    members. This bug affected computations of irradiances and ratios so
    that names used in lists of wavebands where not reflected in the
    output, and the label from the waveband definitions could not be
    overridden as intended by design.

-   Ensure that `irrad()` and `response()` always label the returned
    values correctly. Values are now always tagged according to the
    units used, even when these are selected by setting R options; i.e.,
    values returned by `irrad()` and `response()` are always labeled
    identically as those returned by `q_irrad()`, `e_irrad()`,
    `q_response()` and `e_response()`.

-   Relax check for w.length range allowing longer wavelengths (IR). In
    addition triggering of an error for unlikely wavelengths has been
    replaced by a warning, emitted only when verbose output is enabled.
    An error is triggered now only for wavelengths \< 1 nm.

-   Update `clean.object_spct()` adding parameter `min.Afr` that can be
    used to set a different target than the default of zero and in
    addition implement `clean.object_mspct()`, which was missing.

    ### Compatibility with 'dplyr' (\>= 1.0.0):

    Some small internal changes were needed to avoid errors in calls to
    'dplyr' methods. From user's perspective as 'dplyr' now seems to
    retain classes derived from tibble, it may be necessary in certain
    cases to disable some checks during calls to dplyr methods on
    spectral objects in users' scripts. We have implemented R option
    `"photobiology.check.spct"` to allow checks to be enabled and
    disabled as needed. This option is also useful for optimizing code
    performance.

------------------------------------------------------------------------

***Code-breaking changes in R for photobiology*** By changing the labels
used for values returned by summary functions, this update can break
older code. Additionally, bug fixes that correct the behaviour of
`irrad()`, `fluence()`, `ratio()` and `response()`, and
`trim_waveband()` can potentially break old code. The order in which
attributes are stored has in some cases changed, and while a fixed order
is not guaranteed or even expected for attributes in the R language,
equality tests can return `FALSE` for objects that are functionally
identical but differ in their structure because of their vintage.

------------------------------------------------------------------------

***Code breaking changes in the tidyverse*** Packages in the tidyverse
are evolving to use package 'vctrs' for their implementation (see
[Tidyverse blog](https://www.tidyverse.org/blog/) for news). 'tibble' (3.0.0) seems to
handle row names differently than previous versions. The 'photobiology'
code is not affected, but user code might be affected. Contrary to
earlier versions, 'tibble' (3.0.0) retains member names in vectors,
which means that `"wl.colors"` column created by `tag()` after the
update to 'tibble' will contain named color values. dplyr (1.0.0), not
yet released, better preserves object attributes, causing some
differences in the returned values, as now columns retain S3 class
attributes. The use of 'dplyr' verb`select()` on spectral objects now
fails, while extraction with `[` and `[[` works as expected. However,
`filter()` and `mutate()` work correctly.

------------------------------------------------------------------------

# photobiology 0.9.30

-   Add function `compare_spct()` for comparisons between pairs of
    spectra based on summaries computed over multiple ranges of
    wavelengths.
-   Add utility method `uncollect()` for extracting all members of a
    collection of spectra.
-   Add utility method `wl_thin()` for reducing the storage size of
    \_spct objects by removing data points in regions with only minor
    features.
-   Add `when_measured()` and `when_measured<-()`, `where_measured()`
    and `where_measured<-()`, `what_measured()` and `what_measured<-()`,
    ho`w_measured()` and `how_measured<-()` as an alternative syntax
    consistent with base R for setting and querying the attributes
    when.measured, where.measured, how.measured and what.measured used
    to store metadata in spectral objects and collections.
-   Add spelling synonyms for all normalization-related methods and
    functions.
-   Add specialization for collections of spectra and spelling synonyms
    for `normalized_difference_ind()`.
-   Revise `summary.generic_spct()` to store the name of the summarized
    spectrum object and revise `print.summary_generic_spct()` to display
    the name of the summarized object.
-   Revise fscale() so that by default it sets the "scaled" attribute
    only when the target value for re-scaling is equal to one.
-   Fix bug in `get_attributes()` methods.
-   Add parameter `"address"` to `setWhereMeasured()`, `revise print()`
    methods for spectra to display address when available.
-   Revise `getWhereMeasured()` to consistently return a data.frame,
    even when a geocode is missing.
-   Add function `na_geocode()`, a constructor for a valid geocode data
    frame with all fields set to `NA`s of correct modes.
-   Revise all logic used for geocodes for consistency in returned and
    set values.

------------------------------------------------------------------------

One of these changes to handling of geocodes is potentially
code-breaking as the returned value for missing geocode data is now a
data.frame containing NAs instead of a bare NA. In contrast, handling of
geocodes stored in objects created with earlier versions of the package
is fully consistent with that for new objects.

------------------------------------------------------------------------

# photobiology 0.9.29

-   Add `clean.object_spct()`, `na.rm.object_spct()` and
    `na.omit.object_spct()` methods, which were missing.
-   Correct small off-range errors in `Ler_leaf.spct` data.

# photobiology 0.9.28

-   Fix bug in `getTimeUnit()`.
-   Add `setHowMeasured()` and `getHowMeasured()` using attribute
    `"how.measured"`.

# photobiology 0.9.27

-   Replace non-ASCII characters in documentation.
-   Improve performance of `sun_angles()` and `day_night()` and
    implement Meuus algorithm for julian day. End result is slightly
    slower performance but much higher precision and broader range of
    dates that are handled correctly.
-   Make sure functions work correctly when a tibble is passed as
    argument to parameter geocode.
-   Values returned by `sun_angles()` and `day_night()` have **changed
    very slightly as a result of the improved julian day calculation.**

# photobiology 0.9.26

-   Move 'tibble' from Imports: to Depends: as newly defined classes are
    derived from those exported by package 'tibble', making visibility
    of private specialized methods always necessary.
-   Revise `cps2irrad()` to ensure that all attributes are copied.
-   Implement `s_mean()`, `s_median()`, `s_range()`, `s_sd()`,
    `s_var()`, `s_sum()`, `s_prod()` and `s_mean_se()` methods for
    collections of spectra that operate across spectra returning
    computed values at each wavelength---i.e. a single computed spectral
    object.
-   Add support for `na.action = "replace"` to methods `na.omit()` and
    `na.exclude()`.
-   Add parameter `na.rm` to `smooth_spct()` methods, with unchanged
    default behaviour.
-   Fix bug in `smooth_spct()`: metadata not copied to returned
    spectrum.
-   Fix bug in `insert_hinges_spct()`: metadata not fully copied to
    returned spectrum.
-   Implement methods `"wexler"` and `"goff.gratch"` for water vapor
    pressure estimates. Add function `water_RH()` to compute relative
    humidity from water VP and air temperature.

# photobiology 0.9.25

-   Fix bug in merging of attributes which was causing errors in
    operations between spectra which had specific differences in their
    attributes.
-   Fix incompatibility with tibble (\>= 2.0.0).

# photobiology 0.9.24

-   **Fix major bug** in internally used method `merge_attributes()`
    which in certain cases triggered fatal errors when using operators
    or other methods taking two spectra as arguments. This also affected
    indirectly package 'ooacquire'.
-   Add method `wls_at_target()` implementing search and interpolation
    of wavelengths at which the spectral quantity matches a target
    value. This method is implemented for all spectral classes and
    collections of spectra, expect `raw_spct`, `chroma_spct` and
    `object_spct` and corresponding collections.
-   Add (experimental) support for spectral measurements where the
    exposure time is dependent on the light source instead of on the
    integration time settings in the spectrometer such as when measuring
    xenon-tube flashes.

# photobiology 0.9.23

-   Bug fix: Ensure that `getWhenMeasured()` always returns an object of
    class `"POSIXct"`.
-   Add `"scale.factor"` as formal parameter to methods `irrad()`,
    `q_irrad()`, `e_irrad()`, `response()`, `q_response()` and
    `e_response()`.
-   Change default value for parameter `idx` in the various summary
    methods for collections of spectra. Also add parameter `idx` to *add
    attributes to tibble* methods and rewrite the *get* methods for
    attributes using a simpler, and possibly faster approach.
-   Add parameter `sep` to function `convolve_each()`.

# photobiology 0.9.22

Remove dependency on 'caTools' which is no longer maintained and is
scheduled to be archived in CRAN.

# photobiology 0.9.21

-   Fix bugs in `copy_attributes()` function (some attributes were not
    copied).
-   Add copying of attributes to conversion functions `cps2irrad()`,
    `cps2Tfr()` and `cps2Rfr()`.
-   Add `merge_attributes()` function and fix bugs in handling of
    metadata attributes in math operators, affecting operations between
    two spectral objects.
-   Improve handling of attributes in objects containing multiple
    spectra in tidy (longitudinal) form. Revise `setMultipleWl()` to
    guess number of copies when `multiple.wl` argument is `NULL`. Add
    `setIdFactor()` and `getIdFactor()`, and add parameter `'idfactor'`
    to `setGenericSpct()` and equivalent functions.
-   Handle gracefully attempts to apply `smooth_spct()` method to
    spectra containing `NA`s. (Warn and return as is.)
-   Add parameter `na.rm` to `peaks()`, `valleys()` and related methods
    and functions.
-   Add parameter `na.rm` to `normalize()` methods.
-   Option `"photobiology.verbose"` is initialized to the value of R's
    option `"verbose"` at the time the package is attached.
-   If `photobiology.verbose == TRUE`, starting from this version the
    presence of `NA`s in spectral data triggers warnings.

# photobiology 0.9.20

-   To improve clarity of user code, define wrappers for methods
    `max()`, `min()`, `range()`, `midpoint()`, and `spread()` with names
    with `wl_` prepended (e.g. `wl_max()`). Rename method `spread()` to
    `expanse()` and deprecate use of `spread()`. In the case of
    `wl_expanse()` and `expanse()` the new names prevent a name clash
    with `tidyr::spread()`.
-   Redefine all '`as.`' class-coercion functions into generic methods,
    and add specializations. Default behaviour remains unchanged.
-   Add class-coercion methods to and from matrix (and transfer
    corresponding functions from package 'photobiologyInOut').
-   Add function `relative_AM()` to compute *relative air mass* from sun
    elevation.
-   Add functions for calculation of saturated water vapour pressure
    from air temperature and dew point from water vapour pressure.

# photobiology 0.9.19

-   Fix bug in `sun_angles()` that resulted in wrong azimuth values
    returned in some circumstances. (Bug introduced in version 0.9.12
    when the astronomy code was overhauled in December 2016.)
-   Fix bug in `day_night()` that resulted in a small time shift when a
    time was passed instead of a date as an argument to parameter
    `date`.
-   Fix bug in `irrad()`, `e_irrad()` and `q_irrad()` methods for
    collections of spectra which was resulting in arguments passed to
    parameter `quantity` being silently ignored.
-   Add explicit support for optional parallel execution of apply
    functions, and methods for collections of spectra. Apparently of no
    use under Windows 10, even for large collections of spectra.
-   Revise documentation.

# photobiology 0.9.18

-   Add calibration_spct class and its methods. Update all operators and
    maths. functions to support this new class. Add corresponding test
    cases.
-   Add calibration_mspct class and its methods.
-   Add functions add_attr2tb(), when_measured2tb(), lon2tb(), lat2tb(),
    geocode2tb() and what_measured2tb().
-   Add a formal parameter "attr2tb" to all summary methods for
    collections of spectra, used to automatically call add_attr2tb() to
    add attributes to their output.
-   Fix bug in calculation of solar time of day which could result in
    values \> 24 h.
-   Fix a bug in the extraction operator for spectral objects (could
    lead to infinite recursion in rare occasions).
-   Fix wrong value returned in some cases by is_tagged() due to a bug
    in the definition of tag().
-   Revise vignettes and documentation.
-   Rebuild some of the example data objects adding metadata attributes.

# photobiology 0.9.17

-   Add function `merge2object_spct()` that constructs an `object_spct`
    from a `filter_spct` and a `reflector_spct`, even if wavelengths
    differ. Add methods `T2Afr()` and `Afr2T()` for converting
    transmittance to absorptance and vice versa.
-   Add method `T2T()` for converting between internal and total
    transmittance when total reflectance is available.
-   Implement maths operators and functions for absorptance (`Afr`).
-   Add convenience functions for setting options for the evaluation of
    an expression: `using_photon()`, `using_energy()`, `using_Tfr()`,
    `using_Afr()`, `using_A()`.
-   Add convenience functions for setting options for maths:
    `photon_as_default()`, `energy_as_default()`, `Tfr_as_default()`,
    `Afr_as_default()`, `A_as_default()`.
-   Add convenience functions for setting options for diagnosis:
    `verbose_as_default()`, `strict_range_as_default()`.
-   Add convenience functions for setting options for computation:
    `wb_trim_as_default()`, `use_cached_mult_as_default()`.
-   Add `unset_user_defaults()` to clear all private options recognized
    by this package.
-   Implement new attribute `"Afr.type"` and methods `setAfrType()` and
    `getAfrType()`. Update `copy_attributes()` function and
    `check_spct()` methods to support this new attribute. (subject to
    changes)
-   Revise `cps2Tfr()` and `cps2Rfr()`: simplify validity tests and fix
    bug in `cps2Rfr()`.
-   Add method `join_mspct()` to convert a collection of spectra into a
    data.frame with multiple columns. This method is mainly intended for
    use when exporting spectral data.
-   Add spectral data for an Arabidopsis leaf as `object_spct`,
    `filter_spct`, and `reflector_spct` objects, and two
    `` raw_mspct` `` objects with the corresponding *raw detector
    counts* from the spectrometer.

# photobiology 0.9.16

Rename S3 method `color()` to `color_of()` to avoid a name clash with an
S4 method. A function `color()` calling method `color_of()` maintains
backwards compatibility, but use of `color()` **is now deprecated**.
Improve performance of `sun_angles()`. Convenience functions returning
individual angles, `sun_azimuth()`, `sun_zenith_angle()`, and
`sun_elevation()`, also benefit. Optimize performance of `day_night()`.
Other functions returning individual times for sun positions,
`sunrise_time()`, `noon_time()`, `sunset_time()`, `day_length()`,
`night_length()`, also benefit. **Performance optimization does not
alter the returned values.**

Revise the documentation corresponding to these functions and the
"astronomy" vignette to better describe how vectorization works for
these functions.

# photobiology 0.9.15

Laxen range test for cps in `check_spct.cps_spct()`, so as to avoid
spurious warnings and/or errors. Make `waveband()` with no arguments
valid. Add `na.omit()` and `na.exclude()` methods for all spectral
classes. Functions `na.pass()` and `na.fail()` as defined in R work as
expected for spectra.

# photobiology 0.9.14

Add methods `isValidIntrDesc()` and `isValidInstrSettings()`. Add
methods `trimIntrDesc()` and `trimInstrSettings()`. Improve handling of
metadata attributes when row-binding spectra with `rbindspct()` so that
all metadata can be restored by `subset2mspct()`. Revise print methods
so that they are compatible with new format of attributes in row-bound
spectra. Fix bugs in print method for summaries of spectra, and make
output match that of the print method for spectra (some metadata was not
printed). Fix bug in `tag.mspct()` method (bad default waveband resulted
in no tagging). Add example data for spectral reflectance. Add test
cases for all new code, and add additional tests cases for old code. Add
example data for white LED bulb as `raw_spct`, `cps_spct` and
`source_spct` objects. Improve documentation.

# photobiology 0.9.13

Add functions `sun_elevation()` and `sun_azimuth()`. Add convenience
functions `trim2overlap()` and `extend2extremes()`. Add `print()`
methods for `instr_desc` and `instr_settings` classes. Change threshold
for automatic use of hinges when calculating summaries, so that hinges
are used by default in more cases. Allow lists of length zero as
argument to w.band, and treat them as equivalent to `NULL`. Make
argument `span = NULL` behave in the same way in peak and valley related
summary functions as it has for some time in statistics in package
'ggspectra'. Fix bug in the range check for counts per second which
triggered warnings too easily. Change scaling of D2 and FEL spectra so
that they are expressed in W m-2 nm-1.

# photobiology 0.9.12 (2016-10-22)

Enhancements: Overhaul the astronomy related functions. Non-backward
compatible changes in their parameters. Code rewritten using Meeus'
astronomical algorithms resulting in better accuracy, much faster
calculation of times for solar noon, sunrise and sunset. In addition the
new algorithm can be used for a broader range of years, limited only by
function julian() in R, which ignores leap seconds. Not fully backward
compatible as returned values differ slightly from those in older
versions, as long as indexing has been done by name. Add classes
"solar_time" and "solar_date", constructor solar_time() and methods
is.solar_time(), is.solar_date(), as.solar_date(), format() and print().
Add conversion function as_tod() returning 'time of day' from a datetime
object. sun_angles() and day_night() now return a data.frame instead of
a list. \* The default for twilight argument is no longer "none", but
instead "sunlight". The string "none" is accepted, and is as before
interpreted as the center of the solar disk at the horizon, without
correction for atmospheric diffraction, which is a solar elevation of
zero degrees. In contrast, "sunlight" is interpreted as the upper rim of
the solar disk at the horizon after approximate correction for
atmospheric diffraction, which is a solar elevation of -0.833 degrees,
resulting in longer days.

Revise cbindspct() so that it stores the name of "idfactor" as an
attribute and set the retrieved value of the attribute as default in
subset2mspct(). This makes the conversion from a long form spectral
object back into a collection of spectra 'automatic' in many additional
use cases.

Fix bugs: normalize() now sets scaled attribute to FALSE. fscale() now
sets normalized attribute to FALSE. Allow all summary calculations on
re-scaled data, with a warning. Fix bug in operator '\*' between
reflector.spct and source.spct, which triggered an error.

# photobiology 0.9.11

Fix bug in print.generic_spct().

# photobiology 0.9.10 (not submitted to CRAN)

Argument range.check = NA and range.check = NULL now trigger a message
instead of skipping the test altogether. Logical values retain earlier
behaviour. Fix bug in set functions for spectral classes: option
"photobiology.strict.range" was not obeyed. Revise irrad(), e_irrad()
and q_irrad() to accept as input scaled or normalized source_spct
objects as long as returned quantity is not "total" or "mean", or its
synonym "average". If the returned value will be in relative units
linear scaling or normalization will not invalidate the result, so
earlier restriction was unnecessarily strict. Fix bug in cps2irrad().

# photobiology 0.9.9 (2016-08-16)

Fix bug in rbindspct() affecting binding of object_spct objects. Fix bug
in sunrise_time() and sunset_time() with vector of dates input at polar
region latitudes. Fix bug in e_ratio(). Convert User Guide from PDF
(Rnw) to HTML (Rmd). Fix some minor mistakes in documentation and build
a documentation web site with package 'staticdocs'.

# photobiology 0.9.8

Fix bug: insert_hinges() would fail for "raw_spct" and "cps_spct"
objects. Fix bug: leap year for year 2000. Test, improve and document
function sun_angles(). Clean code for handling of "range" arguments.
Make more consistent handling of "w.band" argument in spectral summary
functions such as irrad(): a numeric range is now a valid argument which
is converted on-the-fly into a waveband object. Add "guard code" to
summary functions triggering errors or warnings when the spct object
passed as argument contains data for multiple spectra in long form. Fix
bug in stepsize.generic_spct() and improve print() output for spct
objects containing data for multiple spectra (in long form) and their
summaries. Improve how print() handles large collections of spectra.

# photobiology 0.9.7 (2016-04-??)

Edit to track changes to the behaviour of functions in package
'lubridate' which triggered errors in several use cases of package
'photobiology'.

# photobiology 0.9.6 (2016-04-02)

Optimize R code for performance: high-level functions faster than in
version 0.9.4 (last version using Rcpp and C++). Add support for
"what.measured" attribute. Allow multiple "countsXXXX" columns in
"raw_spct" objects. Allow multiple "cpsXXXX" columns in "cps_spct"
objects. Fix bug in getInstrSettings(). Allow small rounding and
instrument errors to pass validity checks. Fix bug in setRawSpct().
Implement clean(), normalize(), fscale() and fshift() methods for
"generic_spct", "raw_spct" and "cps_spct" objects. 

New functions: cps2irrad(), cps2Tfr() and cps2Rfr(). Preliminary versions.
Constructors for \_spct objects gain a ... formal argument which allows
addition of arbitrary columns to the objects created. Reorganize
documentation into fewer help files. Fix bug leading to loss of special
attributes. Add test cases. Fix bug in "extract" operator [.object_spct.
Update for compatibility with dplyr (\>= 0.4.3.9001). Update tests for
testthat (\>= 0.11.0.9000).

# photobiology 0.9.5 (2016-02-03)

Remove buggy C++ code and replace it with R code. Performance cost: 20
to 30% slower in high-level functions.

# photobiology 0.9.4

Prepare for CRAN submission. Rename check() -\> check_spct() to avoid
name clash with 'devtools'. Fix bug in clean() methods for collections
of spectra.

# photobiology 0.9.3

Cosmetic changes.

# photobiology 0.9.2

Change tag.generic_spct() and wb2rect_spct() so that "wb.color" is also
added as a column to make it easier to define new ggplot2 stats based on
them. Make tagging compatible with raw_spct and cps_spct objects. Change
behaviour of w.band = NULL as argument to using a single waveband
covering the whole wavelength range of the spectral object, for
consistency with other methods. In addition change behaviour so that the
same columns are added independently of the value supplied as argument
for "w.band", filling the unused columns with NAs when w.band = NA.

# photobiology 0.9.1

Rename some "low level" functions and their formal parameters.

Rename the default methods for midpoint() and stepsize() as methods for
numeric and add a new default method that handles unsupported classes
more gracefully.

Revise color() method to work correctly in borderline cases. Returned
value is always of class character (character() or NA_character\_ in
some cases). The default for the formal parameter 'type' is now
consistently "CMF", instead of in a few cases "both".

Revise interpolate_spct() (interpolate_wl.generic_spct()) to work
correctly in borderline cases and to handle gracefully numeric(). Also
returned value always has the same class and columns, irrespective of
the values of other arguments.

# photobiology 0.9.0

This version removes all two-way dependencies between this package and
other packages in the r4photobiology suite. This should make
installation of the source distribution a lot easier.

## Fix minor bugs

All spectral object constructors return an object with length equal to
zero when called without arguments, or with argument "w.length" of
length equal to zero. Earlier behaviour was to return NA, which is not
consistent with R's usual expectations.

All constructors for collections of spectra return an object with length
equal to zero when called without arguments or with a list of length
zero as first argument "l". Earlier behaviour was to trigger an error,
which is not consistent with R's usual expectations.

Function trim_waveband() was not returning the values described in some
cases. Now the function really behaves as originally documented.

Several functions and methods have a 'range' argument, but semantics was
not identical in all. Now they all use the same semantics, similar to
that used in recent versions of ggplot2 for xlim() and ylim(). Setting a
boundary to NA means that the range extends to include all available
data on a given 'end'. Whenever it makes sense a NULL range means that
the range includes the whole available range of wavelengths. The default
value of 'range' for the recently added method 'fshift()' has changed to
be the first 10 nm at the short wavelength end of the spectrum.

## Removed functionality

Function calc_filter_multipliers() has been deleted.

## New data

Add example data to make this package's documentation independent of
other packages in the suite.

Filters: clear.spct, opaque.spct, polyester.spct, yellow.gel.spct
Sensors: photodiode.spct, ccd.spct Objects: black_body.spct,
white_body.spct, clear_body.spct

## New functionality

Add method clip_wl() to more easily select data for a range of
wavelengths. Add method trim_wl() from trimming, and method
interpolate_wl(). All three methods are implemented for individual
spectra and for collections of spectra. Methods clip_wl() and trim_wl()
have been also implemented for individual waveband objects and for lists
of waveband objects.

## New functionality (only useful to package developers)

Add new classes raw_spct and raw_mspct to store raw counts from
spectrometers, constructors raw_spct() and raw_mspct(), and
corresponding setRawSpct() and is.raw_spct() and as.raw_spct()
functions.

Add methods getInstrDesc() and getInstrSettings(), and the corresponding
setInstrDesc() and setInstrSettings() to handle metadata about
instrument and instrument settings used to acquire the spectral data.
These classes and methods are building blocks to be used in the
development of spectral-data acquisition packages.

## Why "0.9.x" started?

Changes are backwards compatible but new classes were added and several
existing methods were implemented for them. Seven new methods were added
and implemented for all spectral classes and collections of spectra
classes. The new classes are being used for writing a new package
intended to replace package MayaCalc.

# photobiology 0.8.11

Add method clean() for removing out-of-range spectral observations. Add
method fshift() for shifting the zero of the scale used to express
spectral data. Add 'geocode' formal parameter to all functions related
to sun angles.

Clean-up code and vignette. Update vignette.

# photobiology 0.8.10

Add support for "when.measured" and "where.measured" attributes in
generic_spct objects. Add methods setWhereMeasured(),
getWhereMeasured(), setWhenMeasured() and getWhenMeasured(). The first
two methods are compatible with the output of ggmap::geocode() and the
last two with POSIXct objects as returned by many functions in package
lubridate, such as lubridate::now() and lubridate::today(). All these
methods are implemented for both generic_spct and generic_mspct classes.

Implement methods to 'get' and 'set' attributes of summary objects of
spectra.

Revise print() method for spectra to include the new attributes when
available. Rewrite summary() methods and print() method for spectral
summaries. Now output parallels the output of the print() method for
spectra.

# photobiology 0.8.9

Small edit to msdply() to allow setting of column names. Update to
\_mspct summary methods to use the names in the list of wavebands, if
present, as column names for the returned data frame.

# photobiology 0.8.8

Add function trim_mspct(). Implement normalize(), fscale(), peaks(),
valleys(), tag() and untag() methods for generic_mspct and derived
classes. Implement dim\<- for generic_mspct and derived classes. Fix
bugs in msmsply() and mslply(). Implement sign(), round(), signif(),
floor(), ceiling() methods for generic_mspct and derived classes.
Implement log2(), cos(), sin(), tan(), acos(), asin(), acos(), cospi(),
sinpi() and tanpi() methods for generic_mspct and derived classes.

Make internal changes needed for an update to photobiologyInOut. Add
option "photobiolgy.strict.range" to allow changing the default
behaviour for out-of-range values. (Used internally to avoid redundant
checks and repeated warning messages.)

Update package Title in docs. Update User Guide, tables and their
formatting.

# photobiology 0.8.7

Many bug fixes, mostly minor. Rename f_mspct() to msdply() to respect
name conventions in common use. Add mslply() and msaply(). Implement
operators %/% and %%, and function abs() for spectra. Update
documentation. Revise titles of help files.

# photobiology 0.8.6

Bug fixes related to imports from other packages. Rename mutate_mspct()
to msmsply() to respect name conventions in common use. Add constructor
chroma_spct(). Update documentation.

# photobiology 0.8.5

Implement extract and replacement methods for collections of spectra.
(Please, let me know if you encounter any errors with these methods, as
some testthat tests give infinite recursion errors that I cannot
reproduce outside testthat.)

Implement combine method c() for collections of spectra.

Add function convolve_each().

Add helper function shared_member_class().

# photobiology 0.8.4

Expand sun.spct and sun.daily.spct down to 280 nm with zeros.

Add methods for construction of collections of spectra from data frames
with different spectra in side-by-side columns ('untidy' or wide data):
split2source_mspct(), split2response_mscpt(), split2filter_mspct(),
split2reflector_mspct(), and split2cps_mspsct()

Add method for construction of collections of spectra from data frames
with spectral data for different spectra in a single column ('tidy' or
long data), and spectral objects containing several spectra, such as
those returned by function rbindspct(): subset2mspct()

Add print() method for collections of spectra.

Fix minor bug in trim_spct() which was inserting hinges on head and/or
tail expansion even if use.hinges was set to FALSE.

Fix bug in rare borderline cases where NA was being returned instead of
empty spectral objects (objects with zero rows).

# photobiology 0.8.3

Improve handling of multiple.wl \> 1 to avoid spurious warnings. Fix
minor bugs in spread().

# photobiology 0.8.2

Fix bugs in trim_spct() and reduce trimming rounding protection in
operators for spectra.

# photobiology 0.8.1

tag() and untag() by default use copy semantics. Spectral objects with
zero rows are handled cleanly. \$\<- defined for spectra. Fix bugs.

# photobiology 0.8.0

No longer use data.table as a base class for spectral objects. \*\*
Given the size of spectral data the advantages were too limited compared
to the complications introduced. \*\*

Extract and replacement methods "[" and "[\<-" and the subset() function
should now work as expected when applied to spectral objects!
Subscripting of spectra can be used without any restrictions.

Argument passing and assignment semantics follows normal R semantics of
passing of making, \`lazily', a copy, except in a few exceptions.

Improved accuracy of returned values (from the earlier level of
approximately 4 significant digits to close to 7 or 8 significant digits
when using hinges) and improved consistency of returned values between
functions and operators.

New print method for spectral objects.

Several rather minor bugs and inconsistencies fixed, specially in
operators.

Several dependencies and suggests removed. No imports are now visible to
users, which should avoid clashes.

\*\* code is not yet optimized for performance \*\*

\*\* Hopefully not many bugs have been introduced \*\*

# photobiology 0.7.1

Fix bug in new code added in 0.7.0. Update rbinspct() to work well with
\_mspct objects as argument. Add is. and as. functions for all \_mspct
classes. Add function msmsply() and improve code of function f_mspct()
using package plyr. Implement e2q() and q2e() methods for objects of
classes source_mspct and response_mspct. Implement T2A() and A2T()
methods for objects of class filter_mspct. Implement min(), max(),
range(), spread(), midpoint() and stepsize() for objects of class
generic_mspct and derived classes. Implement dim() for objects of class
generic_mspct and derived classes, which also makes available functions
ncol() and nrow() for objects of these classes. Add function
rmDerivedMspct(). Change returned value of rmDerivedSpct() into a
character vector containing the removed class attributes. Update and
reorganize User Guide. Add tests for new code.

# photobiology 0.7.0

Make ratio functions into methods, and add 'wb.trim' as a parameter.
Remove 'pc.out' parameter from methods that had it.

Add \_mspct classes and summary methods for collections of spectral
objects (the 'm' in \_mspct comes from 'multi'. The idea for these new
classes came from a question from Susan Holmes after my talk at
UseR!2015 in Aalborg. She asked whether my package could handle
hyperspectral image data. I answer not, but that this feature could be
easily added. As with the rest of the package I have now a first design
that avoids as much as possible making assumptions about the data, and
allows as much extensibility as possible.

Add attribute 'spct.version' to all spectral objects, and function
getObjectVersion() to query it.

Changed upgrade() method into upgrade_spct() function, and now it also
adds the attribute 'spct.version' during object update.

Add function normalized_diff_ind() usable for calculation of indexes
similar to NDVI (normalized difference vegetation index) for reflectance
and any other spectral summary quantity and two waveband objects. This
was inspired by a poster in Aalborg about package hsdar for analysis of
hyperspectral images.

Revise User Guide and Upgrade Guide.

# photobiology 0.6.8

Functions sunrise_time(), noon_time(), sunset_time(), day_length(),
night_length() and day_night() gain formal parameter 'unit.out', but its
default value preserves earlier behaviour.

New generic functions peaks() and valleys().

subset() methods for spectral objects gain a new parameter 'idx'. During
sub-setting, tags are removed from the returned spectrum.

Functions is.photon_based(), is.energy_based(), is.absorbance_based()
and is.transmittance_based() have been renamed by replacing the dot with
an underscore for consistency with other functions not test for object
class.

Fixed bug in trim_waveband().

# photobiology 0.6.7

Fix bug in functions day_length() and night_length() added in 0.6.6. The
bug was triggering an error in calculations for polar regions. (Test
cases added.)

# photobiology 0.6.6

Package 'lubridate' is no longer imported, so to use functions from this
package, users must explicitly load it.

A user reported that day_night() was not vectorized as stated in the
documentation. This bug has been fixed, and now day_night() and the new
functions noon_time(), sunrise_time(), sunset_time(), day_length() and
night_length() are all vectorized for the 'date' parameter. The default
'tz' for all the functions in now "UTC", to match the default time zone
used to obtain the default for 'date'. The 't' argument of day_night()
has been renamed 'date' for clarity.

These changes, may break old user code! (but needed changes will be
minor).

# photobiology 0.6.5

Add support for "twilight" definition as angle in degrees in function
day_night(), and different angles for sunrise and sunset.

Fix bug in calc_multipliers() that resulted in NAs in cases where
normalization wavelengths fell outside the range of a waveband object
but could be anyway computed without errors because the BSWF was defined
over a wider range of wavelengths than the waveband itself.

Updated dependency versions.

# photobiology 0.6.4

Add support for lubridate::duration: Parameter "time.unit" in methods
and functions related for source_spct and response_spect objects accepts
durations.

New function convertTimeUnit() can be used to modify the time.unit
attribute, re-expressing the spectral data using the new "time.unit".

Old function setTimeUnit() issues a warning when it is used to override
an already set "time.unit" attribute.

Methods irrad(), e_irrad() and q_irrad(), and response(), e_response()
and q_response() gain parameter 'time.unit' which can be used to obtain
their result expressed on a different time-unit basis.

New methods fluence(), e_fluence() and q_fluence().

Fixed bug in absorptance.filter_spct().

Replaced "." by "\_" in the names of the classes returned by summary()
methods for spectra.

# photobiology 0.6.3

Fixed long-standing bug in twilight calculation in function day_night().
irrad(), e_irrad(), and q_irrad() modified so that NAs in the spectral
data outside the range of waveband(s) do not cause the result to be an
NA.

Documentation and User Guide updated.

# photobiology 0.6.2

Added new class cps_spct for uncalibrated spectral data in counts per
second. These are not raw instrument counts, but data linerized and
re-expressed per unit time in seconds. Maths operators and functions are
supported.

Added new accepted value for attribute time.unit: "exposure", meaning
the total exposure time, in which case instead of spectral irradiances,
the objects store spectral fluence or spectral dose.

Added data for CIE standard A illuminant.

Changed package data files to contain each only one R object.

Updated User Guide adding summary tables for spectral objects and their
variables and attributes.

# photobiology 0.6.1

Added User Guide vignette with upgrade instructions.

Added new functions is.old_spct() and upgrade_spectra()

# photobiology 0.6.0

I am reading the new book "R Packages", authored by Hadley Wickham, and
I have found several things to improve in the package, especially with
naming conventions and documentation. I have earlier followed "common
practice" or tried to make names consistent with data.table. As the
package has not yet been publicly released, I am taking this last chance
I have of getting the naming conventions improved.

Major backwards-incompatible changes to function names to avoid method
dispatch problems. All function and class names use underscores. Only S3
methods have dots, and only between method name and class names.
Although many names have changed, they should be still easily
recognizable to current users.

Added test for non-unique wavelength values during spectral object
construction, and added a new formal argument to check() and set...()
functions to allow setting the maximum number of copies allowed. Current
version only triggers a warning.

Generic function Rescale() has been replaced by methods for the generic
function scale() defined in R base, and is_rescaled() has been renamed
is_scaled().

IMPORTANT: because the .spct classes have been renamed, objects created
with earlier versions of the package are not recognized. They need to be
"upgraded" with function upgrade().

Full overhaul of the documentation, using features of latest version of
package roxygen2. Documentation for related functions is now on the same
file, and most pages are cross- referenced.

# photobiology 0.5.19

Added specialized subset() methods for spectral objects as
subset.data.table drops attributes. The new methods are wrappers that
copy the attributes used by this package after the subset() operation is
done by subscripting, including comment, but still drop any attributes
added by users.

Edited check.generic.spct() to generate a warning instead of an error
for spectra of length zero, as can be returned by subset when the
'subset' condition is FALSE for all rows.

# photobiology 0.5.18

Fixed bug in rbindspct().

# photobiology 0.5.17-1

Updated version requirements for all dependencies.

# photobiology 0.5.17

New options "photobiology.use.hinges" (NULL, TRUE, FALSE) with default
NULL and "photobiology.auto.hinges.limit" (wavelength step in nm) with
new default 0.5 nm. Default behaviour is changed both for a bug in
irrad() functions where the limit was 1.1 nm instead of 0.7 nm and as a
result of changing to a smaller default value. Existing and new options
"photobiology.waveband.trim", "photobiology.use.cached.mult" are now
used throughout.

Bug fix: in version 0.5.? the behaviour of the product of a BSWF and a
source spectrum was changed to return a response spectrum. This was not
consistent with the idea that a BSWF is used for quantifying radiation
rather predicting a biological response. The result should remain as
radiation and continue to be expressed in energy or photon based units.
This has been fixed. Starting from the current version, if the waveband
is used with the '\*' operator is a BSWF (tested with function
is_effective()) the resulting spectrum is tagged as being "effective".
The irrad() methods for source.spct have been updated to 'recognize'
effective spectral irradiance when supplied as input.

rbindspct() modified so that if at least one of the spectra in the list
is of effective irradiance, but not all spectra have been calculated
with the same BSWF and normalization, a factor named BSWF is added with
the BSWFs retrieved from the spectra used as levels. The optional factor
selected through parameter idfactor has now its levels always reflecting
the order of the spectra in the input list, even when using a named
input list of spectra. Several rather minor bugs were fixed, including
improved handling of comments and attributes.

summary() methods for spectra have been updated to report the BSWF that
have been used if source spectra contain effective spectral irradiances.

Long-standing unsolved problem partly fixed: interpolate_spct() now in
most cases applies smoothing to the spectrum before attempting
interpolation if the vector of output wavelengths has length \> 1, but
is sparser than then original spectral data.

A few inconsistencies in the formal parameters among similar functions
were fixed by adding the missing parameters. Also some problems in the
documentation for a few functions were fixed.

User guide updated.

# photobiology 0.5.16

Fixed handling of attributes in insert_spct_hinges() and added support
for object.spct to this function. Also fixed small bugs here and there
in the handling of attributes, to make behaviour more consistent. Added
tests to avoid already set attributes to be overwritten by defaults when
using set...() and as....() functions on existing .spct objects.

Added warnings to default methods of generic functions to easy the
diagnostic of problems.

# photobiology 0.5.15

New functions normalize() and rescale() for spectra. Both functions set
object attributes to flag the spectra that have been modified, and no
longer are expressed in absolute units. In addition, tests were added to
summary functions to disallow use of rescaled or normalized spectral
data as input, with the exception of function integrate_spct().

Function rbindspct() was revised to issue a warning when only some of
the spectral objects are rescaled or normalized.

Changes in the code for handling 'time.unit', 'Tfr.type' and 'Rfr.type'
attributes. Added new functions getTimeUnit(), getTfrType() and
getRfrType().

New test cases were added.

# photobiology 0.5.14

New class and methods added: object.spct, with the corresponding
as.object.spct(), is.object.spct(), setObjetcSpct(), and
check.object.spct() functions. Also functions reflectance() and
transmittance() are implemented. New function absorptance() is
implemented only for object.spct objects. Objects of this class can be
used to store corresponding spectral transmittance and spectral
reflectance values. No operators are defined for this class as they
would be ambiguous. The class attribute needs to be changed to either
filter.spct or reflector.spct before using operators, however, no data
is lost in the process, or written except for the class attribute, so
the class can be changed back to object.spct if needed.

New function added: merge.generic.spct() is a wrapper on
merge.data.table() that sets the correct class to the returned merged
spectra, preserves attributes used with .spct objects and by default
merges by w.length.

# photobiology 0.5.13

New function smooth_spct() for spectra added. Based on the smoothing
code in package MayaCalc and optionally behaving as a wrapper to other
smoothers available in R. Implemented for source.spct, filter.spct,
reflector.spct and response.spct objects. Its interface may change. Now
package caTools is needed.

Operations between .spct objects and numeric vectors, possibly of length
one, now preserve other variables (e.g. ID factors) contained in the
spectral objects.

Fixed bug in reflectance() function.

# photobiology 0.5.12

Bad input error reporting improved for off-range Tfr and Rfr values in
check(). Warning given by T2A() when Tfr = 0 result in A = Inf. Added
na.rm = TRUE when min() and max() are called on spectral data.

The formal arguments added to function rbindspct() in version 0.5.7 have
been changed to more closely follow the development version 1.9.5 of
package data.table. This will break user code that uses the previous
syntax as added a few weeks ago to this package.

Fixed bug in handling of 'quantity' arguments "contribution.pc" and
"relative.pc" in function absorbance() which was affecting plot
annotations.

# photobiology 0.5.11

Changed code of irrad() and set\_\_\_Spct() functions to be able to
handle locked data.table objects such as .SD when using by within
`$$ $$` on spct objects. irrad() copies the spectrum only if needed, and the
sorting key is set to "w.length" only if not already set to this same
value.

# photobiology 0.5.10

Added automatic testing based on package testthat.

Test cases for operators and math functions. Test cases for constructors
of spectral objects.

# photobiology 0.5.9

Updated operators' code so that unary '-' and unary '+' now work for
spectra.

Rewrote the code for math functions, removing redundancy, and adding
support for chroma.spct objects.

Fixed a major bug in binary operators code for non-transitive
operations, triggered when the first operand is numeric.

Added a table to the User Guide summarizing the valid and invalid binary
operations, and a second table listing unary operators and math
functions.

Edited the constructors of spct objects so that they recognize
additional variable names for the input data, and rename them.

Added warning messages informing of forbidden operations between
dissimilar spectral objects and between spectral objects and wavebands
are attempted.

# photobiology 0.5.8

Bug in q2e.response.spct() fixed.

# photobiology 0.5.7

rbindspct() now will add, if requested, a factor 'spct.idx' to the
returned spct object, with a different level for each spectrum in the
input list. rbindspct() now works as expected even if spectra of the
same class contain data stored using different types of quantities, e.g.
energy or photon based, or transmittance vs. absorbance. User Guide text
was edited. Code formatting in chunks improved.

# photobiology 0.5.6

Added function clear_photobio.cache() and renamed the cache itself to
"protect"" it from accidental deletion by making it invisible. Optimized
colour calculations for improved speed and debugged the changes. Revised
User Guide. Fixed bug that was causing response() not to be exported.
Fixed bug in operations between spectra and wavebands and added also a
small safety margin to protect from rounding errors. A few minor bugs
also fixed. Fixed an important bug causing error in plot(sun.spct \*
CIE()) (a difficult one to track).

# photobiology 0.5.5

Changed the default for the trimming of wavebands to TRUE, and added a
global option so that the default can be changed globally by setting the
option: e.g. using options(photobiology.waveband.trim = FALSE)

Also added waveband trimming handling to function tag().

# photobiology 0.5.4

Fixed bug in split_bands() function. Handling of names has changed to
some extent, and handling of list objects supplied as x argument is now
done by recursion, adding some additional flexibility.

The User Guide was updated to include examples of the use of
split_bands().

# photobiology 0.5.3

Added function is_tagged() to easily test if an spectrum has been tagged
for plotting. Added function untag() to delete tag data from spectra.

# photobiology 0.5.2

Fixed bug in operators that caused the 'time.unit' attribute not to
propagate to returned values. Fixed bug in operations between wavebands
and spectra which resulted in the returned spectrum being expanded to
the range of the waveband even if the operand spectrum did not overlap
this range.

The option "photobiology.base.unit" has been renamed
"photobiology.radiation.unit", for consistency with package
photobiologygg and for clarity.

# photobiology 0.5.1

Added missing response.spct() function. Fixed two bugs in new_waveband:
1) wb.label was not set when wb.name was auto-generated! 2) Both hinges
were on the outside of the waveband, giving trouble with wavebands
sharing their limits.

Added new function trim_waveband(). Edited irrad(), response(),
transmittance(), absorbance(), reflectance(), ratio() functions for
spectra to use this function internally, and gaining the wb.trim
parameter, which defaults to FALSE. Fixed small numeric errors in the
creation of the waveband for the whole spectrum.

# photobiology 0.5.0

The default data in spectral objects has changed in the present and the
previous versions. This is to avoid, when possible, conversions when
creating spectral objects. Instead, the input is checked in the
functions accepting spectral objects as input. However, values are never
stored as percentages, only as fractions, as this does not affect later
computations, and when a percent is needed the result of the computation
is just multiplied by 100.

I added absorbance related functions, similar to the transmittance and
reflectance ones: absorbance_spct() and absorbance(). All these
functions, as well as those for irradiances and ratios, and responses
now return values with the "radiation.unit" attribute set. In addition
several edits were done to simplify the code when possible. Irradiance,
response, transmittance, reflectance and absorbance functions have
gained a new parameter 'quantity', that can take one value out of
"total", "contribution" (to total as a fraction), "percent"
(contribution as percent), or "average" (the total divided by the
'spread' or width of the waveband in nm). New function response().

Options for operators and functions added:
options(photobiology.base.unit = "energy")
options(photobiology.base.unit = "photon") These options determine
whether operations are carried using energy based or a photon based
quantities. This affects operations involving source.spct and
response.spct objects when using operators, and the defaults for
"unit.out" for functions accepting spectral objects as arguments.

options(photobiology.filter.qty = "transmittance")
options(photobiology.filter.qty = "absorbance") determine whether
operations are carried using transmittance or absorbance quantities
(this affects operations involving filter.spct objects).

The defaults remain the same as in earlier versions, except for irrad()
which instead of having no default for "unit.out" now uses the one set
by the option, if the option is set. Options can be unset with NULL.
options(photobiology.base.unit = NULL) options(photobiology.filter.qty =
NULL)

Moved functions find_peaks(), get_peaks() and get_valleys() from
photobiologygg to this package.

User Guide thoroughly updated to reflect currently available
functionality, including the options described above.

Fixed bug in trim_spct() which led to wrongly copied attribute for
filter.spct objects.

# photobiology 0.4.10

The code for the operators has been rewritten almost from scratch. The
code for the special functions has been edited. In addition to cleaner
and more compact code, now the order of the operators does no longer
affect the dispatch of the operators for spectral objects. There are
some changes in which combinations of arguments are accepted. The
checking is now more strict with respect to what may make sense from the
point of view of the physics behind the calculations. Failed checks
return NA.

All calls to invisible(), except in the \texttt{set\_\_\_Spct}
functions, have been replaced with return() as the lack of visible
output was confusing. It is now documented in the User Guide how to use
options() to control the length of the print output.

# photobiology 0.4.9

Generic functions e_response() and q_response() added, and debugged.
Changed check() for source.spct and response.spct so that it does not
add "energy" or "photon" based spectral data on building objects, as
this adds unneeded computation and complicates the logic by changing the
data on copy operations (which is surprising and can lead to
difficult-to-trace bugs.) This change can be NOT BACKWARDS compatible in
some cases of user code, but should be fine for other functions defined
in the photobiology packages, as they have been edited to maintain the
same visible behaviour.

Fixed a spurious warning by trim_spct() due to minute differences in
wavelengths.

# photobiology 0.4.8

Tweaked handling of w.high by new_waveband(). Fixed bug in function
waveband() which caused it not to properly create weighted wavebands.
Improved handling of names by split_wavebands() and improved the
documentation.

# photobiology 0.4.7

Added function is_waveband().

# photobiology 0.4.6

Fixed various small bugs. Removed rbindlist(), for spectra rbindspct()
should be used instead.

# photobiology 0.4.5

Added waveband() constructor function. Fixed bug in reflectance_spct()
and reflectance() functions. Fixed bug in class_spct(). Revised the user
guide to reflect all changes to the code done in the last few versions.

# photobiology 0.4.4

Added the remaining as.xxxx.spct() functions. Added object creation
functions source.spct(), filter.spct() and reflector.spct(). Cleaned
code to avoid some warnings.

# photobiology 0.4.3

Fixed bug in functions transmittance_spct() and reflectance_spct(), and
added aliases transmittance() and reflectance().

# photobiology 0.4.2

New function rmDerivedSpct() now used internally for setting xxx.spct
classes consistently by removing derived class attributes when class is
set to parent class.

New as.generic.spct() and as.private.spct() functions for spct objects
to easy use of operators.

Bug in default arguments to rbindspct() fixed.

# photobiology 0.4.01

Fixed a troublesome bug caused by a couple of missing exports. MayaCalc
was broken because of this.

Moved functions D2_spectrum() and FEL_spectrum() to this package from
photobiologyLamps

# photobiology 0.4.00

Started new development branch. Added is.xxx.spct() functions for all
classes defined and is_any_spct() and class_spct() functions to query
class of spectral objects.

# photobiology 0.3.16

Edited wb2xx functions. Added other required variables and set their
values to zero, so that the generic.spct objects returned by the
waveband to spectra functions can be safely converted to any of the
specialized spct spectral objects. This makes then labelling of all
types of spectral plots equally easy.

# photobiology 0.3.15

Never released, changes reversed.

# photobiology 0.3.14

The required component of response.spct objects was renamed from
'response' to 's.e.response'. 'response' is also accepted and renamed
into 's.e.response'. If 's.q.response' is present it is accepted and
maintained. New methods e2q() and q2e() for response.spct class.

# photobiology 0.3.13

Replaced several 'return' statements with 'invisible'. Replaced all
remaining data.frame() calls with data.table().

# photobiology 0.3.12

Changed the name of two of the variables added by wb2rect_spct() to more
meaningful ones.

# photobiology 0.3.11

Added setGenericSpct() as a synonym for setGenSpct() to keep naming
consistent. Added functions wb2tagged_spct() and wb2rect_spct() for
creating a tagged generic.spct objects suitable from a list of
wavebands.

# photobiology 0.3.10

Updated vignette to use new package photobiologyWavebands. Added a
generic function tag() and methods for generic.spct and source.spct
objects. tag() can be used to 'tag' the rows based on a list of
wavebands, or simply with rgb color equivalents. A color() method for
source.spct() was added, which returns the overall colour of the
spectrum as a whole. Fixed a 'bug' in tag() that increased execution
speed a lot. Added CIE2008 human luminous efficiency function data.
Added functions irrad(), q_irrad() and e_irrad() as synonyms for
irrad_spct(), q_irrad_spct() and e_irrad_spct(). Added functions
irrad(), q_ratio(), e_ratio(), eq_ratio(), and qe_ratio() as synonyms
for q_ratio_spct(), e_ratio_spct(), eq_ratio_spct() and qe_ratio_spct().
eq_ratio_spct() is also a new function. Added function split_bands() for
creating lists of unweighted wavebands. Several bug fixes, and new
function put_hinges(). Added function rbindspct() and redefined
rbindlist() (hiding data.table::rbindlist() by exporting
photobiology::rbindlist() instead.) Also some changes to color() with
new methods for numeric and lists.

# photobiology 0.3.9

Added summary() and corresponding print() methods for spectral objects,
and stepsize() default and "generic.spct" methods. Changed
SetSourceSpct() to accept a "time.unit" argument, which is used to set
an attribute in "source.spct" objects (defaults to "second" making it
transparently backwards compatible).

Performance tests show that with default settings the irradiance
functions for "source.spct" objects execute significantly faster than
the vector based ones. Operators for spectra are at the moment
relatively very slow.

# photobiology 0.3.8

Now chromaticity and colour matching function data are as chroma.spct
objects. This required also edits to many of the colour related
functions. Added CMF for honeybees and example daily solar spectral
data. Edited irrad_spct() so that results always include waveband names.

# photobiology 0.3.7

Added new class response.spct and revised the operators' definitions.

# photobiology 0.3.6

The set...Spct() functions until 0.3.5 had different semantics than
setDT in data.table. This has now been fixed. Data for CIE D65
illuminant has been added. The imports and re-exports from data.table
now include additional functions including some new to data.table 1.9.3.
Added new class chroma.spct for storing x y z chromaticy data and the
corresponding functions setChromaSpct and operators \* and /.

# photobiology 0.3.5

Fixed "bug"" in trim_spct(), which was returning data.table objects
instead of objects of the same type as supplied as argument. A bug in
version 1.9.2 of data.table was behind some other problems, so added
requirement for data.table 1.9.3+ to DESCRIPTION.

# photobiology 0.3.4

Added a class for reflectance spectra and updated operator definitions.
Added function setReflectorSpct(). New functions interpolate_spct(),
reflectance_spct(), transmittance_spct(), q_ratio_spct(),
e_ratio_spct(), qe_ratio_spct(), rgb_spct() and check(),
check.generic.spct(), check.filter.spct(), check.reflector.spct(), and
check.source.spct(). Operator \*.generic.spct now accepts wavebands as
second argument. Debugged min(), max(), range(), midpoint(), spread(),
labels() for spectra. Updated setGenericSpct(), setFilterSpct() and
setSourceSpct() to use check() to check and fix if possible missing
variable, or add the missing variables as a vector of NAs.

Updated the vignette with examples of the use of spct objects and the
functions and operators defined for them. Modified packages so as to
consistently use .spct as the name tag of all objects representing
spectra. Recently .dt was used and earlier, .data, that is still a
fall-back in the current version.

TODO: make sure all spct functions and operators work correctly when the
argument spct objects are incomplete.

# photobiology 0.3.3

New functions insert_spct_hinges(), trim_spct(), integrate_spct(),
average_spct(), irrad_spct(), e_irrad_spct(), q_irrad_scpt(),
transmittance_spct(), reflectance_spct(), ... , is_effective.waveband(),
Small edits to new_waveband() Example solar spectrum data now available
as "source.spct" objects and as "data.table" and "data.frame" so as not
to break any old code.

# photobiology 0.3.2

Function day_night had bugs, now fixed.

# photobiology 0.3.1

Added functions sun_angles and day_night.

# photobiology 0.3.0

**MAJOR UPDATE** 

There are some _backward incompatibilities_

Saved wavebands will give an error, but any code used to create them
should still work unchanged.

If you used the parameter use.cpp.code in your code, you will just need
to delete the corresponding arguments. C++ code is always used.

I hope I haven't broken anything else but I have made quite major
internal changes.

At the moment the code is a bit slower than before. I will most likely
be able to cure this at a later stage.

***

The package now depends on data.table which is now imported. As MayaCalc
already returns data.table objects, hopefully this will help optimize
performance and memory use. Users should be aware that data.tables are
copied by reference, and a true copy is needed, function copy should be
used.

I am making the code more object (S3) oriented. Three new classes are
exported "filter.spct", "source.spct", and "generic.spct", and three new
functions were added to set the class of objects to these classes. If
the object is not already a data.table, it is setDT. This allows me to
create a more user-friendly interface by overloading operators. Now the
four basic operations work on spectra, and one can also use these same
operations and \^between a spectrum a numerical vector of length one or
longer. Many maths functions and some summary functions have specialized
versions for spectra. This is my attempt to make calculations intuitive,
but do not expect them to be fast.

Added generic functions A2T and T2A for converting absorbance to
transmittance and back. The generic function takes a numerical argument,
and there are special versions for spectra. When the argument is a
"filter.spc" object they return a copy of the original object with the
calculated quantity added. If the requested output already exists, it is
NOT recalculated, just a copy of the argument is returned, if the input
data is missing, the corresponding output is set to NA.

Added a label field to waveband object to store an optional label for
plotting, which by default is set equal to the name of the waveband.
labels.waveband was changed to return in field label this label instead
of the name, and the name is returned in the field name. Also added a
new function normalization.waveband that returns the normalization
wavelength.

calc_filter_multipliers now accepts both a filter name as a character
string or an object of class "filter.spct" as argument to parameter
filter (which was renamed from filter.name). In addition the function
acquired a pc.in logical parameter, an the pc parameter was renamed
pc.out. This allows any combination of input and output as percentages
or fractions.

# photobiology 0.2.25

Modified calc_filter_multipliers to work with photobiologyFilters 0.1.8
and later Both calc_filter_multipliers and calc_source_multipliers now
find by name both objects with ending .dt and .data. If both exist .dt
is used.

# photobiology 0.2.24

Added functions oper_spectra, prod_spectra, div_spectra, and
subt_spectra, and redefined sum_spectra based on oper_spectra. Behaviour
of sum_spectra has not changed.

# photobiology 0.2.23

Removed center_wl center_wl.generic and center_wl.waveband Added spread
spread.generic and spread.waveband Added midpoint midpoint.generic and
midpoint.waveband

# photobiology 0.2.22

Moved from photobiologyLamps and renamed function calc_lamp_output() to
calc_source_output(). Moved from photobiologyLamps
calc_filter_multipliers(). Added two examples of how to draw Maxwell's
triangles.

# photobiology 0.2.21

Fixed a bug in s.e.irrad2rgb() caused by not explicitly handling the
case when the wavelength range was completely outside outside the range
of the the CMF or CC data. In other words the function crashed when the
colour of invisible radiation was requested. Now "black" is returned in
this case.

# photobiology 0.2.20

New function w_length_range2rgb(), which takes a range of wavelengths,
and returns the equivalent RGB colour assuming equal spectral (energy)
irradiance at all wave lengths within the range. Modified the
calculation of color for wavebands, using this new function.

# photobiology 0.2.19

Added waveband versions of generic functions, min(), max(), range() and
labels(). Updated trim_tails() so that it can also expand a spectrum.
New functions: w_length2rgb() to calculate RGB colour definitions for
monochromatic light. s.e.irrad2rgb() to calculate RGB colour definition
from spectral irradiance data or reflectance data. Data on x, y, z
coordinates for CIE chromaticity coordinates and colour matching
functions added.

# photobiology 0.2.18

Fixed a bug in split_irradiance() and improved handling of cut point
values outside the data range, unsorted cut-points, and redundant
cut-points. For 'good' data there are no changes in output, except in
the case when three cut points were supplied, as the third one was being
ignored. The handling of out-of-range cut points has been changed so
that the endpoints are moved to the extreme wavelengths of the data.

# photobiology 0.2.17

The cache is now created in the "emptyenv". This solves a long-standing
error message during installation.

# photobiology 0.2.16

Changed irradiance() so that it optionally accepts a list of wave_bands
and returns a vector or irradiances. The returned values have a names
attribute set to the names of the elements in the wave_band list. If the
names are missing, it returns the "name"s stored in the wave_band
objects.

wave_band objects are now S3 objects with class "waveband", a
print.waveband() function was added, which is called whenever the
generic print function is called on a "waveband" object.

New convenience functions split_irradiance(), split_energy_irradiance()
and split_photon_irradiance() for integrating unweighted irradiances for
a series of contiguous wavelength ranges.

# photobiology 0.2.15

Updated NAMESPACE to fix error caused by changes in Rcpp or R.

# photobiology 0.2.14 (2014-01-05)

Updated User Guide. Modified interpolate_spectrum so that it uses cubic
spline interpolation when the number of data points is small, and linear
interpolation otherwise.

# photobiology 0.2.13

Moved peaks-related functions to new package photobiologygg and removed
dependence on splus2R.

# photobiology 0.2.12

Added new functions find_peaks and get_peaks, which are wrappers built
on top of function peaks from package splus2R, which becomes a new
dependency.

# photobiology 0.2.11

Added Rcpp version of insert_hinges() using the much faster binary
search algorithm.

# photobiology 0.2.10

Small, backward compatible change to trim_tails: added fill parameter.
If NULL (default) the function works as earlier deleting the tails. If
another value is passed as fill, the s.irrad values in the tails are
replaced with this value.

# photobiology 0.2.9

Added two new functions: interpolate_spectrum() and sum_spectra()

# photobiology 0.2.7 + 0.2.8

The vignette "Examples" was moved to a separate package because it used
a very large data set. An "benign"" error in new_waveband() was fixed
(this has no effect on calculations). The data set sun.data was updated
to contain instantaneous data, and its documentation updated to better
describe the data.

# photobiology 0.2.6 (2013-08-27)

A vignette "Examples" was added. Other examples will be added later on.
No changes to R or C++ code.

# photobiology 0.2.5 (2013-08-11)

A vignette "Manual" was added. No changes to R or C++ code.

# photobiology 0.2.4

As profiling showed that calls to insert_hinges() and within this
function to which() were using most of run time of the functions, a new
parameter was added to functions so that the use of hinges can be
switched on and off. The default is to use hinges for interpolation only
if the wavelength resolution of the spectrum is worse than 1.0 nm. This
causes only small errors in the results of calculations (\< 1%) but
reduces run time very significantly. The defaults should work just fine
in most cases, in which case there is no need to do any changes to user
code. In addition new_waveband() was modified so that if no SWF is used,
then hinges are not included in the waveband defined.

# photobiology 0.2.2-0.2.3

Revised photon_ratio() to use same optimizations.

Improved implementation of cache. Now it is enclosed in an environment
that is created when the package is loaded and removed when it is
unloaded.

# photobiology 0.2.1

Factored to a separate function the code to check spectral data, and
added argument to irradiance to switch off the checking.

Added a "name" field to wavebands and a wb.name argument to
new_waveband(), with a default so that no old code will be broken.

Added argument to calc_multipliers() to enable caching. Disabled by
default.

Added argument to irradiance() to enable/disable checking of spectral
data. Enabled by default.

Added argument to irradiance() to enable use of speed optimizations.

Defaults to arguments ensure that nothing changes from earlier versions
unless new features are explicitly enabled.

# photobiology 0.2.0

Reimplemented irradiance(), calc_multipliers, e2qmol_multipliers(), and
added e2quantum_multipliers(). Also reimplemented new_waveband(),
changing its interface. Some changes in other photobiologyxxxx will be
needed.

# photobiology 0.1.0 (2013-07-08)

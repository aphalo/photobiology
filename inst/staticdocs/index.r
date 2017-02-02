sd_section("Package overview", "",
           c(
             "photobiology-package"
           )
)

sd_section("Constructors",
  "Spectral data is stored in objects of different classes depending on the physical quantities in the data.",
  c(
    "source_spct",
    "as.generic_spct",
    "setGenericSpct",
    "generic_mspct",
    "as.generic_mspct",
    "waveband",
    "split_bands",
    "wb2rect_spct",
    "wb2spct",
    "wb2tagged_spct",
    "solar_time",
    "as.solar_date")
)

sd_section("Class", "",
           c(
             "is.generic_spct",
             "is.generic_mspct",
             "is.waveband",
             "class_spct",
             "is.summary_generic_spct",
             "spct_classes",
             "mspct_classes",
             "summary_spct_classes",
             "rmDerivedSpct",
             "rmDerivedMspct",
             "shared_member_class",
             "is.solar_time",
             "is.solar_date"
           )
)

sd_section("Metadata attributes", "",
           c(
             "is_photon_based",
             "getBSWFUsed",
             "setBSWFUsed",
             "is_effective",
             "getInstrDesc",
             "setInstrDesc",
             "trimInstrDesc",
             "isValidInstrDesc",
             "getInstrSettings",
             "setInstrSettings",
             "trimInstrSettings",
             "isValidInstrSettings",
             "getMultipleWl",
             "setMultipleWl",
             "getNormalized",
             "setNormalized",
             "is_normalized",
             "getRfrType",
             "setRfrType",
             "getScaled",
             "setScaled",
             "is_scaled",
             "getTfrType",
             "setTfrType",
             "is_absorbance_based",
             "getTimeUnit",
             "setTimeUnit",
             "getWhatMeasured",
             "setWhatMeasured",
             "getWhenMeasured",
             "setWhenMeasured",
             "getWhereMeasured",
             "setWhereMeasured",
             "is_tagged",
             "dim.generic_mspct",
             "labels"
           )
)

sd_section("Summaries",
  "It is often desired to calculate different summaries from spectral data.",
  c(
    "summary",
    "irrad",
    "e_irrad",
    "q_irrad",
    "fluence",
    "e_fluence",
    "q_fluence",
    "response",
    "e_response",
    "q_response",
    "absorbance",
    "absorptance",
    "transmittance",
    "reflectance",
    "average_spct",
    "integrate_spct",
    "color",
    "rgb_spct",
    "e_ratio",
    "q_ratio",
    "eq_ratio",
    "qe_ratio",
    "normalized_diff_ind",
    "min",
    "max",
    "range",
    "spread",
    "midpoint",
    "stepsize"
  )
)

sd_section("Spectral features", "",
           c(
             "peaks",
             "valleys"
           )
)

sd_section("Transformations", "",
           c(
             "normalize",
             "normalization",
             "fscale",
             "fshift",
             "clean",
             "insert_hinges",
             "interpolate_wl",
             "interpolate_spct",
             "smooth_spct",
             "trim_wl",
             "trim_spct",
             "clip_wl",
             "convertTimeUnit",
             "e2q",
             "q2e",
             "T2A",
             "A2T",
             "tag",
             "untag",
             "cps2irrad"
           )
)

sd_section("Maths", "",
  c(
    "plus-.generic_spct",
    "div-.generic_spct",
    "minus-.generic_spct",
    "mod-.generic_spct",
    "slash-.generic_spct",
    "times-.generic_spct",
    "^.generic_spct",
    "Trig",
    "log",
    "MathFun",
    "sign",
    "round"
    )
)

sd_section("Extract, replace and combine", "",
           c(
             "Extract",
             "Extract_mspct",
             "c.generic_mspct",
             "rbindspct",
             "merge.generic_spct",
             "split2mspct",
             "subset2mspct"
           )
)

sd_section("Apply", "",
           c(
             "msmsply",
             "convolve_each"
           )
)

sd_section("Astronomy", "",
           c(
             "day_night",
             "sun_angles",
             "solar_time",
             "as_tod",
             "tz_time_diff"
           )
)

sd_section("Output", "",
           c(
             "print",
             "print.summary_generic_spct",
             "print.waveband",
             "print.solar_time",
             "format.solar_time"
           )
)

sd_section("Data", "radiation sources",
           c(
             "sun.spct",
             "sun.data",
             "sun.daily.spct",
             "sun.daily.data",
             "white_led.raw.spct",
             "white_led.cps.spct",
             "white_led.source.spct",
             "A.illuminant.spct",
             "D65.illuminant.spct",
             "D2_spectrum",
             "D2.UV586",
             "D2.UV653",
             "D2.UV654",
             "FEL_spectrum",
             "FEL.BN.9101.165"
           )
)

sd_section("Data", "objects including filters",
           c(
             "yellow_gel.spct",
             "polyester.spct",
             "clear.spct",
             "opaque.spct",
             "clear_body.spct",
             "white_body.spct",
             "black_body.spct"
           )
)

sd_section("Data", "light sensors",
           c(
             "ccd.spct",
             "photodiode.spct"
           )
)

sd_section("Data", "human and animal vision",
           c(
             "ciev10.spct",
             "ciev2.spct",
             "ciexyzCC10.spct",
             "ciexyzCC2.spct",
             "ciexyzCMF10.spct",
             "ciexyzCMF2.spct",
             "beesxyzCMF.spct"
           )
)

sd_section("Housekeeping", "",
           c(
             "getSpctVersion",
             "getMspctVersion",
             "is.old_spct",
             "check_spct",
             "check_w.length",
             "checkTimeUnit",
             "copy_attributes",
             "normalize_range_arg",
             "upgrade_spct",
             "upgrade_spectra"
           )
)

sd_section("Low-level functions for numeric vectors (deprecated for user code)",
           "These functions are for use by developers when maximum flexibility is needed. They are not meant for everyday use.",
           c(
             "as_energy",
             "as_quantum",
             "as_quantum_mol",
             "calc_multipliers",
             "calc_source_output",
             "check_spectrum",
             "div_spectra",
             "e2qmol_multipliers",
             "e2quantum_multipliers",
             "e2quantum_multipliers",
             "energy_irradiance",
             "energy_ratio",
             "find_peaks",
             "get_peaks",
             "insert_spct_hinges",
             "integrate_xy",
             "interpolate_spectrum",
             "irradiance",
             "oper_spectra",
             "photon_irradiance",
             "photon_ratio",
             "photons_energy_ratio",
             "prod_spectra",
             "s_e_irrad2rgb",
             "split_energy_irradiance",
             "split_irradiance",
             "split_photon_irradiance",
             "subt_spectra",
             "sum_spectra",
             "trim_tails",
             "trim_waveband",
             "w_length_range2rgb",
             "w_length2rgb",
             "waveband_ratio"
           )
)


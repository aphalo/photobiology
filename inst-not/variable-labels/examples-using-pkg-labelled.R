library(photobiology)
library(labelled)
# current RStudio does not display labels for classes derived from data.frame
# or tibble, only for these exact classes
labelled_sun.spct <- sun.spct
var_label(labelled_sun.spct) <- make_var_labels(labelled_sun.spct)
str(labelled_sun.spct)
labels(labelled_sun.spct)
str(var_label(labelled_sun.spct))
generate_dictionary(labelled_sun.spct)
generate_dictionary(sun.spct)
str(remove_var_label(labelled_sun.spct))

norm_sun.spct <- normalize(sun.spct)
var_label(norm_sun.spct) <- make_var_labels(norm_sun.spct)
str(norm_sun.spct)
labels(norm_sun.spct)
str(var_label(norm_sun.spct))
generate_dictionary(norm_sun.spct)

# If supported we will need to update labels when they become invalid
# variable labels could be supported as optiomal

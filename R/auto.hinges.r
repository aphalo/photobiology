auto_hinges <- function(w.length,
                        step.limit = getOption("photobiology.auto.hinges.limit", default = 0.5)) {
  use.hinges <- stepsize(w.length)[2] > step.limit
}

auto_hinges <- function(w.length,
                        step.limit = getOption("photobiology.auto.hinges.limit", default = 0.5)) {
  #                        step.limit = 0.5) {
  stepsize(w.length)[2] > step.limit
}

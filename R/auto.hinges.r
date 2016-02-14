auto_hinges <- function(w.length,
                        step.limit = 0.5) {
  stepsize(w.length)[2] > step.limit
}

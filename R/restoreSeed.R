# Restores the global RNG state saved before a function reseeded it.
# oldSeed is the previous value of .Random.seed, or NULL if the RNG was not yet initialised.
restoreSeed <- function(oldSeed) {
  if (is.null(oldSeed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", oldSeed, envir = globalenv())
  }
}

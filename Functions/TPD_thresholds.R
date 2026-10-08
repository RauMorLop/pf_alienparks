
  # This function creates the isolines that will enclosure density clusters
  # in the TPD function

tpd_thresholds <- function(z, probs = c(.95, .50, .20, .05)) {
  z <- z[is.finite(z)]
  total <- sum(z)
  z <- sort(z / total, decreasing = TRUE)
  cum <- cumsum(z)
  out <- vapply(
    probs,
    function(p) z[min(which(cum >= p))] * total,
    numeric(1)
  )
  names(out) <- c("t95", "t50", "t20", "t05")
  out
}
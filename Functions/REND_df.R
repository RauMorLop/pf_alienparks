  # This function transforms the REND output into a dataframe.
  # It allows us to manage Functional diversity metrics and calculate
  # their SES easier

REND_df <- function(x) {
  
  level <- intersect(c("group", "species"), names(x))[1]
  
  out <- x[[level]]
  
  metric_names <- c(
    FRichness = "Richness",
    FEvenness = "Evenness",
    FDivergence = "Divergence",
    Richness = "Richness",
    Evenness = "Evenness",
    Divergence = "Divergence"
  )
  
  metrics <- intersect(names(out),names(metric_names))
  
  dplyr::bind_rows(
    lapply(metrics,
      function(metric) {
        data.frame(
          Group = names(out[[metric]]),
          Metric = unname(metric_names[metric]),
          Obs = as.numeric(out[[metric]])
        )
      }
    )
  )
}
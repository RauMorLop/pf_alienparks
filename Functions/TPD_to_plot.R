
  # This function calls for data in the TPD object that can be used in ggplot()

tpd_to_plot <- function(x) {
  grid <- as.data.frame(x$data$evaluation_grid)
  names(grid)[1:2] <- c("Axis1", "Axis2")
  dplyr::bind_rows(
    lapply(
      names(x$TPDs),
      function(group) {
        data.frame(
          Group = group,
          grid,
          Density = x$TPDs[[group]]
        )
      }
    )
  )
}
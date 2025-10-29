#' @importFrom magrittr %>%
NULL

#' Create timeline visualization
#'
#' @param timeline_obj MSKTimeline object
#' @param patients Character vector of patient IDs to plot (max 20)
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' plot_timeline(timeline_data)
#' }
plot_timeline <- function(timeline_obj, patients = NULL) {
  
  check_required_packages("ggplot2", "plotting")
  
  if (!methods::is(timeline_obj, "MSKTimeline")) {
    stop("timeline_obj must be an MSKTimeline object", call. = FALSE)
  }
  
  # Combine all timeline data with proper column alignment
  all_data <- list()
  all_colnames <- list()
  
  for (name in names(timeline_obj@data)) {
    df <- timeline_obj@data[[name]]
    df$timeline_source <- name
    all_data[[name]] <- df
    all_colnames[[name]] <- names(df)
  }
  
  if (length(all_data) == 0) {
    stop("No data available for plotting", call. = FALSE)
  }
  
  # Get union of all column names
  all_cols <- unique(unlist(all_colnames))
  
  # Align all data frames to have the same columns
  aligned_data <- lapply(all_data, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    
    # Add missing columns with NA
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    
    # Reorder columns to match
    df[, all_cols, drop = FALSE]
  })
  
  # Combine
  plot_data <- do.call(rbind, aligned_data)
  rownames(plot_data) <- NULL
  
  # Check for required columns
  if (!"DAYS" %in% names(plot_data)) {
    stop("Timeline data must have DAYS column for plotting", call. = FALSE)
  }
  
  if (!"PATIENT_ID" %in% names(plot_data)) {
    stop("Timeline data must have PATIENT_ID column for plotting", call. = FALSE)
  }
  
  # Remove rows with missing DAYS
  plot_data <- plot_data[!is.na(plot_data$DAYS), ]
  
  if (nrow(plot_data) == 0) {
    stop("No data available for plotting after removing NA values", call. = FALSE)
  }
  
  # Filter patients if specified
  if (!is.null(patients)) {
    plot_data <- plot_data[plot_data$PATIENT_ID %in% patients[1:min(20, length(patients))], ]
  } else {
    # Take first 10 patients
    unique_patients <- unique(plot_data$PATIENT_ID)[1:min(10, length(unique(plot_data$PATIENT_ID)))]
    plot_data <- plot_data[plot_data$PATIENT_ID %in% unique_patients, ]
  }
  
  if (nrow(plot_data) == 0) {
    stop("No data available for selected patients", call. = FALSE)
  }
  
  # Use EVENT_TYPE for coloring if available, otherwise use timeline_source
  color_var <- if ("EVENT_TYPE" %in% names(plot_data)) "EVENT_TYPE" else "timeline_source"
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = DAYS, y = PATIENT_ID)) +
    ggplot2::geom_point(ggplot2::aes_string(color = color_var), size = 3, alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Patient Timeline Visualization",
      subtitle = paste("Showing", length(unique(plot_data$PATIENT_ID)), "patients"),
      x = "Days from Baseline",
      y = "Patient ID",
      color = gsub("_", " ", tools::toTitleCase(color_var))
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 8)
    )
  
  return(p)
}
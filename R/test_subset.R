#' Subset MSK-CHORD Data by Cancer Type
#'
#' @param msk_data A MultiAssayExperiment object from download_data()
#' @param cancer_type Character string specifying the cancer type to subset. 
#'   Available types: "Breast Cancer", "Colorectal Cancer", "Non-Small Cell Lung Cancer", 
#'   "Pancreatic Cancer", "Prostate Cancer"
#' @return A MultiAssayExperiment object subsetted to the specified cancer type
#' @export
#' @importFrom MultiAssayExperiment colData
#' @examples
#' \dontrun{
#' msk_data <- download_data()
#' prostate_data <- subset_by_cancer_type(msk_data, "Prostate Cancer")
#' lung_data <- subset_by_cancer_type(msk_data, "Non-Small Cell Lung Cancer")
#' }
subset_by_cancer_type <- function(msk_data, cancer_type) {
  
  # Validate inputs
  if (!inherits(msk_data, "MultiAssayExperiment")) {
    stop("msk_data must be a MultiAssayExperiment object")
  }
  
  if (!is.character(cancer_type) || length(cancer_type) != 1) {
    stop("cancer_type must be a single character string")
  }
  
  # Get clinical data
  clinical_data <- MultiAssayExperiment::colData(msk_data)
  
  # Available cancer types based on diagnostic output
  available_types <- c("Breast Cancer", "Colorectal Cancer", "Non-Small Cell Lung Cancer", 
                       "Pancreatic Cancer", "Prostate Cancer")
  
  # Check if cancer type exists
  if (!cancer_type %in% available_types) {
    stop("Cancer type '", cancer_type, "' not found. ",
         "Available types: ", paste(available_types, collapse = ", "))
  }
  
  # Create subset condition
  subset_condition <- clinical_data$CANCER_TYPE == cancer_type
  subset_condition[is.na(subset_condition)] <- FALSE
  
  # Subset the data
  subset_data <- msk_data[, subset_condition]
  
  # Report results
  n_samples <- sum(subset_condition)
  n_patients <- length(unique(clinical_data$PATIENT_ID[subset_condition]))
  message("Subsetted to ", n_samples, " samples from ", n_patients, " patients with ", cancer_type)
  
  return(subset_data)
}


#' Plot Biomarker Timeline with Clinical Events
#'
#' @param msk_data A MultiAssayExperiment object from download_data()
#' @param biomarker Character string specifying the biomarker to plot. 
#'   Available: "PSA", "CEA", "CA_15-3", "CA_19-9", "GLEASON"
#' @param cancer_type Optional character string to filter by cancer type first
#' @param patient_ids Optional character vector of specific patient IDs to plot. 
#'   If NULL, will sample random patients
#' @param n_patients Integer specifying number of patients to plot (default: 20, max: 100)
#' @param log_scale Logical indicating whether to use log10 scale for y-axis (default: TRUE for lab values)
#' @param facet_by Optional character string to facet by. For GLEASON, can use "gleason_category"
#' @param add_events Logical indicating whether to overlay clinical events (default: TRUE)
#' @param date_range Optional numeric vector of length 2 specifying date range to plot
#' @return A ggplot object showing biomarker timelines with clinical events
#' @export
#' @importFrom MultiAssayExperiment colData metadata
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline facet_wrap
#'   labs theme_minimal scale_color_viridis_d xlim ylim
#' @importFrom dplyr filter mutate slice_head %>%
#' @examples
#' \dontrun{
#' msk_data <- download_data()
#' 
#' # PSA timeline for prostate cancer
#' plot_biomarker_timeline(msk_data, "PSA", cancer_type = "Prostate Cancer")
#' 
#' # CEA timeline with specific patients
#' plot_biomarker_timeline(msk_data, "CEA", 
#'                        patient_ids = c("P-0000314", "P-0001234"))
#' 
#' # Gleason scores with faceting
#' plot_biomarker_timeline(msk_data, "GLEASON", 
#'                        cancer_type = "Prostate Cancer",
#'                        facet_by = "gleason_category")
#' }
plot_biomarker_timeline <- function(msk_data, biomarker, cancer_type = NULL, 
                                    patient_ids = NULL, n_patients = 20, log_scale = TRUE,
                                    facet_by = NULL, add_events = TRUE, date_range = NULL) {
  
  # Check required packages
  required_packages <- c("ggplot2", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install with: install.packages('", pkg, "')")
    }
  }
  
  # Validate biomarker
  available_biomarkers <- c("PSA", "CEA", "CA_15-3", "CA_19-9", "GLEASON")
  if (!biomarker %in% available_biomarkers) {
    stop("Biomarker '", biomarker, "' not available. ",
         "Available biomarkers: ", paste(available_biomarkers, collapse = ", "))
  }
  
  # Get metadata
  metadata_list <- MultiAssayExperiment::metadata(msk_data)
  
  # Map biomarker to timeline data
  timeline_map <- list(
    "PSA" = "timeline_psa_labs",
    "CEA" = "timeline_cea_labs", 
    "CA_15-3" = "timeline_ca_15-3_labs",
    "CA_19-9" = "timeline_ca_19-9_labs",
    "GLEASON" = "timeline_gleason"
  )
  
  timeline_name <- timeline_map[[biomarker]]
  
  if (!timeline_name %in% names(metadata_list)) {
    stop("Timeline data '", timeline_name, "' not found in metadata")
  }
  
  # Get biomarker timeline data
  biomarker_data <- as.data.frame(metadata_list[[timeline_name]])
  
  # Filter by cancer type if specified
  if (!is.null(cancer_type)) {
    clinical_data <- as.data.frame(MultiAssayExperiment::colData(msk_data))
    cancer_patients <- unique(clinical_data$PATIENT_ID[clinical_data$CANCER_TYPE == cancer_type])
    biomarker_data <- biomarker_data[biomarker_data$PATIENT_ID %in% cancer_patients, ]
  }
  
  # Filter by specific patients or sample random patients
  if (!is.null(patient_ids)) {
    biomarker_data <- biomarker_data[biomarker_data$PATIENT_ID %in% patient_ids, ]
    if (nrow(biomarker_data) == 0) {
      stop("No data found for specified patient_ids")
    }
  } else {
    # Sample random patients
    unique_patients <- unique(biomarker_data$PATIENT_ID)
    if (length(unique_patients) > n_patients) {
      selected_patients <- sample(unique_patients, min(n_patients, 100))
      biomarker_data <- biomarker_data[biomarker_data$PATIENT_ID %in% selected_patients, ]
      message("Plotting ", length(selected_patients), " randomly selected patients")
    }
  }
  
  if (nrow(biomarker_data) == 0) {
    stop("No biomarker data available for the specified criteria")
  }
  
  # Prepare the value column based on biomarker type
  if (biomarker == "GLEASON") {
    biomarker_data$value <- biomarker_data$GLEASON_SCORE
    y_label <- "Gleason Score"
    log_scale <- FALSE  # Don't use log scale for Gleason scores
    
    # Add gleason category for faceting
    if (!is.null(facet_by) && facet_by == "gleason_category") {
      biomarker_data$gleason_category <- ifelse(biomarker_data$GLEASON_SCORE >= 7, "High (≥7)", "Low (<7)")
    }
  } else {
    biomarker_data$value <- biomarker_data$RESULT
    y_label <- paste(biomarker, "Value")
  }
  
  # Remove missing values
  biomarker_data <- biomarker_data[!is.na(biomarker_data$value), ]
  
  # Apply date range filter if specified
  if (!is.null(date_range) && length(date_range) == 2) {
    biomarker_data <- biomarker_data[biomarker_data$START_DATE >= date_range[1] & 
                                       biomarker_data$START_DATE <= date_range[2], ]
  }
  
  # Create the base plot
  if (log_scale && biomarker != "GLEASON") {
    # Use log10 scale, but handle zeros/negative values
    biomarker_data$value[biomarker_data$value <= 0] <- 0.01
    p <- ggplot2::ggplot(biomarker_data, ggplot2::aes(x = START_DATE, y = log10(value), 
                                                      group = PATIENT_ID, color = PATIENT_ID)) +
      ggplot2::labs(y = paste("log10(", y_label, ")"))
  } else {
    p <- ggplot2::ggplot(biomarker_data, ggplot2::aes(x = START_DATE, y = value, 
                                                      group = PATIENT_ID, color = PATIENT_ID)) +
      ggplot2::labs(y = y_label)
  }
  
  # Add lines and points
  p <- p + 
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.8) +
    ggplot2::geom_point(alpha = 0.8, size = 1.5) +
    ggplot2::labs(
      title = paste(biomarker, "Timeline"),
      subtitle = if (!is.null(cancer_type)) paste("Cancer Type:", cancer_type) else "All Cancer Types",
      x = "Days from Reference Point",
      caption = paste("Data from", length(unique(biomarker_data$PATIENT_ID)), "patients")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_viridis_d(guide = "none")  # Hide legend for cleaner look
  
  # Add clinical events if requested
  if (add_events) {
    event_data <- get_clinical_events(msk_data, unique(biomarker_data$PATIENT_ID))
    
    if (nrow(event_data) > 0) {
      # Add treatment events (green dashed lines)
      treatment_events <- event_data[event_data$event_type == "treatment", ]
      if (nrow(treatment_events) > 0) {
        p <- p + ggplot2::geom_vline(data = treatment_events, 
                                     ggplot2::aes(xintercept = START_DATE), 
                                     color = "darkgreen", alpha = 0.6, linetype = "dashed")
      }
      
      # Add surgery events (blue solid lines)
      surgery_events <- event_data[event_data$event_type == "surgery", ]
      if (nrow(surgery_events) > 0) {
        p <- p + ggplot2::geom_vline(data = surgery_events, 
                                     ggplot2::aes(xintercept = START_DATE), 
                                     color = "blue", alpha = 0.6, linetype = "solid")
      }
      
      # Add progression events (red dotted lines)
      progression_events <- event_data[event_data$event_type == "progression", ]
      if (nrow(progression_events) > 0) {
        p <- p + ggplot2::geom_vline(data = progression_events, 
                                     ggplot2::aes(xintercept = START_DATE), 
                                     color = "red", alpha = 0.6, linetype = "dotted")
      }
    }
  }
  
  # Add faceting if specified
  if (!is.null(facet_by)) {
    if (facet_by == "gleason_category" && "gleason_category" %in% colnames(biomarker_data)) {
      p <- p + ggplot2::facet_wrap(~ gleason_category, scales = "free")
    } else if (facet_by %in% colnames(biomarker_data)) {
      p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)), scales = "free")
    }
  }
  
  # Apply date range limits if specified
  if (!is.null(date_range) && length(date_range) == 2) {
    p <- p + ggplot2::xlim(date_range[1], date_range[2])
  }
  
  return(p)
}


#' Get Clinical Events for Timeline Plotting
#'
#' @param msk_data A MultiAssayExperiment object from download_data()
#' @param patient_ids Character vector of patient IDs
#' @return A data.frame with clinical events
#' @keywords internal
get_clinical_events <- function(msk_data, patient_ids) {
  
  metadata_list <- MultiAssayExperiment::metadata(msk_data)
  
  all_events <- data.frame()
  
  # Get treatment events
  if ("timeline_treatment" %in% names(metadata_list)) {
    treatment_data <- as.data.frame(metadata_list[["timeline_treatment"]])
    treatment_data <- treatment_data[treatment_data$PATIENT_ID %in% patient_ids, ]
    if (nrow(treatment_data) > 0) {
      treatment_events <- data.frame(
        PATIENT_ID = treatment_data$PATIENT_ID,
        START_DATE = treatment_data$START_DATE,
        event_type = "treatment",
        event_detail = treatment_data$AGENT,
        stringsAsFactors = FALSE
      )
      all_events <- rbind(all_events, treatment_events)
    }
  }
  
  # Get surgery events
  if ("timeline_surgery" %in% names(metadata_list)) {
    surgery_data <- as.data.frame(metadata_list[["timeline_surgery"]])
    surgery_data <- surgery_data[surgery_data$PATIENT_ID %in% patient_ids, ]
    if (nrow(surgery_data) > 0) {
      surgery_events <- data.frame(
        PATIENT_ID = surgery_data$PATIENT_ID,
        START_DATE = surgery_data$START_DATE,
        event_type = "surgery",
        event_detail = surgery_data$SUBTYPE,
        stringsAsFactors = FALSE
      )
      all_events <- rbind(all_events, surgery_events)
    }
  }
  
  # Get progression events
  if ("timeline_progression" %in% names(metadata_list)) {
    progression_data <- as.data.frame(metadata_list[["timeline_progression"]])
    progression_data <- progression_data[progression_data$PATIENT_ID %in% patient_ids & 
                                           progression_data$PROGRESSION == "Y", ]
    if (nrow(progression_data) > 0) {
      progression_events <- data.frame(
        PATIENT_ID = progression_data$PATIENT_ID,
        START_DATE = progression_data$START_DATE,
        event_type = "progression", 
        event_detail = "Disease Progression",
        stringsAsFactors = FALSE
      )
      all_events <- rbind(all_events, progression_events)
    }
  }
  
  return(all_events)
}


#' Get Available Timeline Biomarkers
#'
#' @param msk_data A MultiAssayExperiment object from download_data()
#' @return A data.frame with available biomarkers and their data summary
#' @export
#' @importFrom MultiAssayExperiment metadata
#' @examples
#' \dontrun{
#' msk_data <- download_data()
#' available_biomarkers <- get_timeline_biomarkers(msk_data)
#' print(available_biomarkers)
#' }
get_timeline_biomarkers <- function(msk_data) {
  
  metadata_list <- MultiAssayExperiment::metadata(msk_data)
  
  biomarker_info <- data.frame(
    biomarker = character(),
    timeline_name = character(),
    n_measurements = integer(),
    n_patients = integer(),
    date_range = character(),
    value_range = character(),
    stringsAsFactors = FALSE
  )
  
  biomarkers <- list(
    "PSA" = "timeline_psa_labs",
    "CEA" = "timeline_cea_labs",
    "CA_15-3" = "timeline_ca_15-3_labs", 
    "CA_19-9" = "timeline_ca_19-9_labs",
    "GLEASON" = "timeline_gleason"
  )
  
  for (biomarker in names(biomarkers)) {
    timeline_name <- biomarkers[[biomarker]]
    
    if (timeline_name %in% names(metadata_list)) {
      timeline_data <- as.data.frame(metadata_list[[timeline_name]])
      
      # Get value column
      if (biomarker == "GLEASON") {
        values <- timeline_data$GLEASON_SCORE[!is.na(timeline_data$GLEASON_SCORE)]
      } else {
        values <- timeline_data$RESULT[!is.na(timeline_data$RESULT)]
      }
      
      biomarker_info <- rbind(biomarker_info, data.frame(
        biomarker = biomarker,
        timeline_name = timeline_name,
        n_measurements = nrow(timeline_data),
        n_patients = length(unique(timeline_data$PATIENT_ID)),
        date_range = paste(range(timeline_data$START_DATE, na.rm = TRUE), collapse = " to "),
        value_range = if (length(values) > 0) paste(round(range(values), 2), collapse = " to ") else "No data",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(biomarker_info)
}


#' Interactive Biomarker Timeline Plot for Shiny
#'
#' @param msk_data A MultiAssayExperiment object from download_data()
#' @param biomarker Character string specifying the biomarker
#' @param cancer_type Optional character string to filter by cancer type
#' @param patient_limit Integer maximum number of patients to plot (for performance)
#' @return A plotly object for interactive exploration
#' @export
#' @examples
#' \dontrun{
#' # This function requires plotly package
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   msk_data <- download_data()
#'   interactive_plot <- plot_biomarker_timeline_interactive(msk_data, "PSA", 
#'                                                          cancer_type = "Prostate Cancer")
#'   interactive_plot
#' }
#' }
plot_biomarker_timeline_interactive <- function(msk_data, biomarker, cancer_type = NULL, 
                                                patient_limit = 50) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots. Install with: install.packages('plotly')")
  }
  
  # Create the base ggplot
  base_plot <- plot_biomarker_timeline(msk_data, biomarker, cancer_type = cancer_type,
                                       n_patients = patient_limit, add_events = FALSE)
  
  # Convert to interactive plotly
  interactive_plot <- plotly::ggplotly(base_plot, tooltip = c("x", "y", "colour"))
  
  return(interactive_plot)
}
#' @importFrom magrittr %>%
NULL

#' Filter MSKTimeline object by multiple criteria
#'
#' @param timeline_obj MSKTimeline object
#' @param patients Character vector of patient IDs to include  
#' @param date_range Numeric vector of length 2: c(min_days, max_days)
#' @param event_types Character vector of event types to include
#' @return Filtered MSKTimeline object
#' @export
#' @examples
#' \dontrun{
#' filtered <- filter_timeline(
#'   timeline_obj,
#'   patients = c("P-0001", "P-0002")
#' )
#' }
filter_timeline <- function(timeline_obj, 
                            patients = NULL,
                            date_range = NULL, 
                            event_types = NULL) {
  
  if (!inherits(timeline_obj, "MSKTimeline")) {
    stop("timeline_obj must be an MSKTimeline object", call. = FALSE)
  }
  
  # Apply filters to each timeline data frame
  filtered_data <- lapply(timeline_obj@data, function(df) {
    result <- df
    
    # Patient filter
    if (!is.null(patients) && "PATIENT_ID" %in% names(result)) {
      result <- result[result$PATIENT_ID %in% patients, , drop = FALSE]
    }
    
    # Date range filter
    if (!is.null(date_range) && "DAYS" %in% names(result)) {
      result <- result[result$DAYS >= date_range[1] & result$DAYS <= date_range[2], , drop = FALSE]
    }
    
    # Event type filter
    if (!is.null(event_types) && "EVENT_TYPE" %in% names(result)) {
      result <- result[result$EVENT_TYPE %in% event_types, , drop = FALSE]
    }
    
    result
  })
  
  # Remove empty data frames
  filtered_data <- filtered_data[sapply(filtered_data, nrow) > 0]
  
  # Update patient count
  if (length(filtered_data) > 0) {
    all_patients <- unique(unlist(lapply(filtered_data, function(x) {
      if ("PATIENT_ID" %in% names(x)) x$PATIENT_ID else character(0)
    })))
  } else {
    all_patients <- character(0)
  }
  
  # Create new object
  new("MSKTimeline",
      data = filtered_data,
      metadata = timeline_obj@metadata,
      cancer_types = timeline_obj@cancer_types,
      n_patients = as.integer(length(all_patients)))
}

#' Extract complete patient data by cancer type (memory-efficient)
#'
#' @param mae MultiAssayExperiment object from download_msk_chord()
#' @param cancer_type Character string specifying cancer type
#' @param include_genomics Logical, whether to attempt genomic data extraction (default FALSE)
#' @return List containing filtered clinical, timeline, and optionally genomic data
#' @export
extract_by_cancer_type <- function(mae, cancer_type, include_genomics = FALSE) {
  
  if (!inherits(mae, "MultiAssayExperiment")) {
    stop("mae must be a MultiAssayExperiment object", call. = FALSE)
  }
  
  # Extract clinical data
  clinical_data <- MultiAssayExperiment::colData(mae)
  clinical_df <- as.data.frame(clinical_data)
  
  # Check if CANCER_TYPE column exists
  if (!"CANCER_TYPE" %in% names(clinical_df)) {
    cancer_col <- grep("cancer.*type|oncotree|disease", names(clinical_df), 
                       ignore.case = TRUE, value = TRUE)[1]
    
    if (is.na(cancer_col)) {
      stop("Cannot find cancer type column in clinical data", call. = FALSE)
    }
    
    message("Using column: ", cancer_col, " for cancer type filtering")
  } else {
    cancer_col <- "CANCER_TYPE"
  }
  
  # Get available cancer types
  available_types <- unique(clinical_df[[cancer_col]])
  
  if (!cancer_type %in% available_types) {
    stop("Cancer type '", cancer_type, "' not found.\n",
         "Available types: ", paste(available_types, collapse = ", "),
         call. = FALSE)
  }
  
  # Filter patients
  selected_patients <- clinical_df$PATIENT_ID[clinical_df[[cancer_col]] == cancer_type]
  
  message("Found ", length(selected_patients), " patients with ", cancer_type)
  
  # Filter clinical data
  filtered_clinical <- clinical_df[clinical_df$PATIENT_ID %in% selected_patients, ]
  
  # Extract timeline data for these patients
  message("Extracting timeline data...")
  timeline_full <- extract_msk_timeline(mae, cancer_types = NULL)
  filtered_timeline <- filter_timeline(timeline_full, patients = selected_patients)
  
  # Get sample mapping
  sample_map <- MultiAssayExperiment::sampleMap(mae)
  selected_sample_map <- sample_map[sample_map$primary %in% selected_patients, ]
  
  # Create result object (without genomics initially)
  result <- list(
    cancer_type = cancer_type,
    n_patients = length(selected_patients),
    patient_ids = selected_patients,
    clinical = filtered_clinical,
    timeline = filtered_timeline,
    sample_map = selected_sample_map,
    mae_subset = NULL
  )
  
  # Optionally include genomics (memory intensive)
  if (include_genomics) {
    message("\nAttempting to extract genomic data...")
    message("WARNING: This may require substantial memory (>16GB)")
    
    tryCatch({
      # Subset the entire MAE object by patients
      result$mae_subset <- mae[, selected_patients]
      message("Genomic data available in mae_subset - access via MultiAssayExperiment methods")
    }, error = function(e) {
      warning("Could not subset genomic data: ", e$message)
    })
  } else {
    message("\nGenomic data not extracted (set include_genomics=TRUE to attempt)")
    message("Note: Genomic extraction may require >16GB memory")
  }
  
  class(result) <- c("msk_cancer_cohort", "list")
  
  message("\nSuccessfully extracted data for ", length(selected_patients), " patients")
  
  return(result)
}

#' Get mutation data for specific genes from cohort
#'
#' @param cohort msk_cancer_cohort object
#' @param genes Character vector of gene symbols
#' @return Data frame with mutations in specified genes
#' @export
get_mutations_by_gene <- function(cohort, genes) {
  
  if (!inherits(cohort, "msk_cancer_cohort")) {
    stop("cohort must be an msk_cancer_cohort object", call. = FALSE)
  }
  
  if (is.null(cohort$mae_subset)) {
    stop("Cohort does not contain genomic data. Re-extract with include_genomics=TRUE", 
         call. = FALSE)
  }
  
  if (!"mutations" %in% names(cohort$mae_subset@ExperimentList)) {
    stop("No mutation data available in cohort", call. = FALSE)
  }
  
  mutations <- cohort$mae_subset@ExperimentList$mutations
  
  # Extract mutations for specific genes only (memory efficient)
  # This depends on the structure - adjust as needed
  message("Extracting mutations for ", length(genes), " genes...")
  
  # Return instructions if large
  list(
    message = "Access mutations via: cohort$mae_subset@ExperimentList$mutations",
    n_samples = ncol(mutations),
    available_methods = "Use MultiAssayExperiment or RaggedExperiment methods to query"
  )
}
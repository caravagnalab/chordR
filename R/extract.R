#' @importFrom magrittr %>%
NULL

#' Extract timeline data from MSK-CHORD dataset
#'
#' @param mae MultiAssayExperiment object from download_msk_chord()
#' @param cancer_types Character vector of cancer types to include (optional)
#' @return MSKTimeline object
#' @export
#' @examples
#' \dontrun{
#' msk_data <- download_msk_chord()
#' timeline_data <- extract_msk_timeline(msk_data)
#' }
extract_msk_timeline <- function(mae, cancer_types = NULL) {
  
  message("Extracting timeline data...")
  
  # Get metadata
  metadata_obj <- MultiAssayExperiment::metadata(mae)
  
  # Extract timeline configurations
  timeline_configs <- get_default_timeline_configs()
  
  # Extract all timelines
  timeline_data <- extract_all_timelines(metadata_obj, timeline_configs)
  
  # Create MSKTimeline object
  result <- create_timeline_object(timeline_data, mae)
  
  message("Successfully extracted timeline data for ", 
          result@n_patients, " patients")
  
  return(result)
}

#' Get default timeline configurations
#' @return List of timeline configurations
#' @keywords internal
get_default_timeline_configs <- function() {
  list(
    labs = list(
      fields = c("timeline_ca_15-3_labs", "timeline_ca_19-9_labs", 
                 "timeline_cea_labs", "timeline_psa_labs"),
      event_type = "LAB_RESULT"
    ),
    treatments = list(
      fields = c("timeline_treatment", "timeline_surgery", "timeline_radiation"),
      event_type = "TREATMENT"
    )
  )
}

#' Extract all timeline types
#' @param metadata_obj Metadata from MultiAssayExperiment
#' @param timeline_configs List of timeline configurations
#' @return Named list of timeline data frames
#' @keywords internal
extract_all_timelines <- function(metadata_obj, timeline_configs) {
  
  timeline_data <- lapply(timeline_configs, function(config) {
    extract_single_timeline(metadata_obj, config)
  })
  
  # Remove NULL entries
  timeline_data[!sapply(timeline_data, is.null)]
}

#' Standardize timeline column names
#' @param data Data frame from timeline
#' @return Data frame with standardized columns
#' @keywords internal
standardize_timeline_columns <- function(data) {
  
  # Rename START_DATE to DAYS if it exists
  if ("START_DATE" %in% names(data) && !"DAYS" %in% names(data)) {
    names(data)[names(data) == "START_DATE"] <- "DAYS"
  }
  
  # Rename STOP_DATE to END_DAYS if needed
  if ("STOP_DATE" %in% names(data) && !"END_DAYS" %in% names(data)) {
    names(data)[names(data) == "STOP_DATE"] <- "END_DAYS"
  }
  
  # Rename TEST to TEST_NAME if needed
  if ("TEST" %in% names(data) && !"TEST_NAME" %in% names(data)) {
    names(data)[names(data) == "TEST"] <- "TEST_NAME"
  }
  
  # Rename RESULT to VALUE if needed
  if ("RESULT" %in% names(data) && !"VALUE" %in% names(data)) {
    names(data)[names(data) == "RESULT"] <- "VALUE"
  }
  
  return(data)
}

#' Extract single timeline type
#' @param metadata_obj Metadata from MultiAssayExperiment
#' @param config Timeline configuration list
#' @return Data frame with timeline data or NULL
#' @keywords internal
extract_single_timeline <- function(metadata_obj, config) {
  
  available_fields <- intersect(config$fields, names(metadata_obj))
  
  if (length(available_fields) == 0) {
    return(NULL)
  }
  
  # Combine data from multiple fields with column alignment
  all_data <- list()
  all_colnames <- list()
  
  for (field in available_fields) {
    data <- metadata_obj[[field]]
    
    if (!is.null(data) && nrow(data) > 0) {
      # Convert to regular data frame
      data <- as.data.frame(data)
      
      # Standardize column names
      data <- standardize_timeline_columns(data)
      
      # Add metadata columns
      data$EVENT_TYPE <- config$event_type
      data$SOURCE_FIELD <- field
      
      all_data[[field]] <- data
      all_colnames[[field]] <- names(data)
    }
  }
  
  if (length(all_data) == 0) {
    return(NULL)
  }
  
  # If only one data frame, return it
  if (length(all_data) == 1) {
    return(all_data[[1]])
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
  
  # Now rbind will work
  combined_data <- do.call(rbind, aligned_data)
  rownames(combined_data) <- NULL
  
  return(combined_data)
}

#' Create timeline object from extracted data
#' @param timeline_data List of timeline data frames
#' @param source_mae Original MultiAssayExperiment object
#' @return MSKTimeline object
#' @keywords internal
create_timeline_object <- function(timeline_data, source_mae) {
  
  # Extract clinical data
  clin_data <- MultiAssayExperiment::colData(source_mae)
  
  # Calculate patient count
  if (length(timeline_data) > 0) {
    all_patients <- unique(unlist(lapply(timeline_data, function(x) {
      if ("PATIENT_ID" %in% names(x)) x$PATIENT_ID else character(0)
    })))
  } else {
    all_patients <- character(0)
  }
  
  # Create metadata
  metadata <- list(
    creation_date = Sys.Date(),
    source_study = attr(source_mae, "study_id") %||% "msk_chord_2024"
  )
  
  # Create object
  new("MSKTimeline",
      data = timeline_data,
      metadata = metadata,
      cancer_types = character(0),  # Would need to extract from clinical data
      n_patients = as.integer(length(all_patients)))
}

#' Prepare complete INCOMMON dataset with calculated depths
#'
#' @param cohort msk_cancer_cohort object
#' @param mae Original MultiAssayExperiment
#' @param max_samples Maximum samples to extract
#' @param min_purity Minimum tumor purity threshold (default 0.1, i.e. 10%)
#' @return Clean dataset ready for INCOMMON
#' @export
prepare_incommon_data <- function(cohort, mae, max_samples = 100, min_purity = 0.1) {
  
  cat("PREPARING INCOMMON DATASET\n")
  cat("==========================\n\n")
  
  # Get mutation samples
  mut_samples <- cohort$sample_map$colname[cohort$sample_map$assay == "mutations"]
  
  if (length(mut_samples) > max_samples) {
    cat("Limiting to first", max_samples, "of", length(mut_samples), "mutation samples\n")
    mut_samples <- mut_samples[1:max_samples]
  }
  
  cat("Extracting mutations from", length(mut_samples), "samples...\n")
  mutations <- mae@ExperimentList$mutations[, mut_samples]
  
  # Convert to data frame
  mut_gr <- as(mutations, "GRangesList")
  
  mutation_list <- list()
  samples_processed <- 0
  
  for (i in seq_along(mut_gr)) {
    sample_id <- names(mut_gr)[i]
    
    if (length(mut_gr[[i]]) > 0) {
      df <- as.data.frame(mut_gr[[i]])
      df$SAMPLE_ID <- sample_id
      
      # extract gene symbols from rownames
      df$Hugo_Symbol <- names(mut_gr[[i]])
      
      # Get patient ID and purity
      patient_id <- cohort$sample_map$primary[cohort$sample_map$colname == sample_id][1]
      df$PATIENT_ID <- patient_id
      
      purity <- cohort$clinical$TUMOR_PURITY[cohort$clinical$PATIENT_ID == patient_id][1]
      purity <- purity /100 # Convert to fraction
      df$PURITY <- purity  

      
      mutation_list[[sample_id]] <- df
      samples_processed <- samples_processed + 1
    }
    
    if (samples_processed %% 10 == 0) {
      cat("  Processed", samples_processed, "samples...\n")
    }
  }
  
  mutations_df <- bind_rows_safe(mutation_list)
  cat("Total mutations extracted:", nrow(mutations_df), "\n\n")
  
  # Calculate depths if missing
  cat("Calculating sequencing depths...\n")
  if (!"t_depth" %in% names(mutations_df) || all(is.na(mutations_df$t_depth))) {
    mutations_df$t_depth <- mutations_df$t_alt_count + mutations_df$t_ref_count
    cat("  Calculated t_depth from alt + ref counts\n")
  }
  
  if (!"n_depth" %in% names(mutations_df) || all(is.na(mutations_df$n_depth))) {
    mutations_df$n_depth <- mutations_df$n_alt_count + mutations_df$n_ref_count
    cat("  Calculated n_depth from alt + ref counts\n")
  }
  
  # Calculate VAF
  mutations_df$VAF <- mutations_df$t_alt_count / mutations_df$t_depth
  
  cat("\nFiltering for data quality...\n")
  initial_count <- nrow(mutations_df)
  
  # Filter for complete cases
  mutations_clean <- mutations_df %>%
    dplyr::filter(
      !is.na(t_alt_count),
      !is.na(t_ref_count),
      !is.na(t_depth),
      !is.na(PURITY),
      PURITY >= min_purity,  # Purity is now a fraction (0-1)
      t_depth > 0
    )

  
  cat("  Removed", initial_count - nrow(mutations_clean), "mutations with missing/invalid data\n")
  cat("  Retained", nrow(mutations_clean), "high-quality mutations\n")
  
  # Extract CNA data
  cat("\nExtracting copy number data...\n")
  cna_samples <- intersect(unique(mutations_clean$SAMPLE_ID),
                           cohort$sample_map$colname[cohort$sample_map$assay == "cna"])
  
  cna_data <- NULL
  if (length(cna_samples) > 0) {
    cna <- mae@ExperimentList$cna
    
    if (inherits(cna, "SummarizedExperiment")) {
      cna_matrix <- SummarizedExperiment::assays(cna)[[1]][, cna_samples]
      cna_data <- as.data.frame(as.matrix(cna_matrix))
      cna_data$GENE <- rownames(cna_matrix)
      cat("  Extracted CNA for", length(cna_samples), "samples\n")
    }
  }
  
  # Create summary statistics
  cat("\n=== FINAL DATASET SUMMARY ===\n")
  cat("Samples:", length(unique(mutations_clean$SAMPLE_ID)), "\n")
  cat("Patients:", length(unique(mutations_clean$PATIENT_ID)), "\n")
  cat("Mutations:", nrow(mutations_clean), "\n")
  cat("Genes mutated:", length(unique(mutations_clean$Hugo_Symbol)), "\n")
  
  cat("\nData completeness:\n")
  cat("  t_alt_count: 100%\n")
  cat("  t_ref_count: 100%\n")
  cat("  t_depth: 100%\n")
  cat("  Purity: 100%\n")
  
  cat("\nData ranges:\n")
  cat("  Purity:", paste(round(range(mutations_clean$PURITY), 3), collapse = "-"), "(fraction)\n")
  cat("  t_depth:", paste(range(mutations_clean$t_depth), collapse = "-"), "\n")
  cat("  VAF:", paste(round(range(mutations_clean$VAF), 3), collapse = "-"), "\n")
  
  if (!is.null(cna_data)) {
    cat("  CNA values:", paste(range(cna_data[,-ncol(cna_data)], na.rm = TRUE), collapse = " to "), "\n")
  }
  
  # Create result object
  result <- list(
    mutations = mutations_clean,
    cna = cna_data,
    purity = unique(mutations_clean[, c("PATIENT_ID", "SAMPLE_ID", "PURITY")]),
    parameters = list(
      n_samples = length(unique(mutations_clean$SAMPLE_ID)),
      n_patients = length(unique(mutations_clean$PATIENT_ID)),
      n_mutations = nrow(mutations_clean),
      min_purity = min_purity
    )
  )
  
  class(result) <- c("incommon_ready", "list")
  return(result)
}

#' Print method for incommon_ready
#' @export
print.incommon_ready <- function(x, ...) {
  cat("INCOMMON-Ready Dataset\n")
  cat("======================\n")
  cat("Samples:", x$parameters$n_samples, "\n")
  cat("Patients:", x$parameters$n_patients, "\n")
  cat("Mutations:", x$parameters$n_mutations, "\n")
  cat("\nReady for INCOMMON Bayesian inference\n")
  invisible(x)
}

#' Filter mutations by gene list
#'
#' @param incommon_data Object from prepare_incommon_data()
#' @param genes Character vector of gene symbols to filter for
#' @return Filtered incommon_ready object
#' @export
#' @examples
#' \dontrun{
#' # Filter for common cancer genes
#' filtered <- filter_by_genes(incommon_data, c("TP53", "KRAS", "EGFR"))
#' }
filter_by_genes <- function(incommon_data, genes) {
  
  if (!inherits(incommon_data, "incommon_ready")) {
    stop("Input must be an incommon_ready object from prepare_incommon_data()")
  }
  
  if (!"Hugo_Symbol" %in% names(incommon_data$mutations)) {
    stop("Hugo_Symbol column not found in mutations data")
  }
  
  cat("Filtering mutations for", length(genes), "genes...\n")
  
  # Filter mutations
  filtered_muts <- incommon_data$mutations %>%
    dplyr::filter(Hugo_Symbol %in% genes)
  
  cat("  Retained", nrow(filtered_muts), "of", nrow(incommon_data$mutations), "mutations\n")
  
  # Get samples that still have mutations
  remaining_samples <- unique(filtered_muts$SAMPLE_ID)
  cat("  Retained", length(remaining_samples), "of", 
      incommon_data$parameters$n_samples, "samples\n")
  
  # Filter CNA data if present
  filtered_cna <- NULL
  if (!is.null(incommon_data$cna)) {
    cna_genes <- intersect(genes, incommon_data$cna$GENE)
    if (length(cna_genes) > 0) {
      filtered_cna <- incommon_data$cna %>%
        dplyr::filter(GENE %in% genes)
      
      # Keep only columns for remaining samples
      sample_cols <- intersect(remaining_samples, names(filtered_cna))
      filtered_cna <- filtered_cna %>%
        dplyr::select(GENE, dplyr::all_of(sample_cols))
      
      cat("  Retained CNA data for", length(cna_genes), "genes\n")
    }
  }
  
  # Filter purity data
  filtered_purity <- incommon_data$purity %>%
    dplyr::filter(SAMPLE_ID %in% remaining_samples)
  
  # Create filtered result
  result <- list(
    mutations = filtered_muts,
    cna = filtered_cna,
    purity = filtered_purity,
    parameters = list(
      n_samples = length(remaining_samples),
      n_patients = length(unique(filtered_muts$PATIENT_ID)),
      n_mutations = nrow(filtered_muts),
      min_purity = incommon_data$parameters$min_purity,
      filtered_genes = genes
    )
  )
  
  class(result) <- c("incommon_ready", "list")
  return(result)
}

#' Verify data completeness for INCOMMON analysis
#'
#' @param incommon_data Object from prepare_incommon_data()
#' @return Invisible TRUE if all checks pass, stops with error otherwise
#' @export
verify_incommon_data <- function(incommon_data) {
  
  cat("VERIFYING INCOMMON DATA COMPLETENESS\n")
  cat("====================================\n\n")
  
  required_mut_cols <- c("Hugo_Symbol", "SAMPLE_ID", "PATIENT_ID", 
                         "t_alt_count", "t_ref_count", "t_depth", 
                         "VAF", "PURITY")
  
  # Check mutations data
  cat("Checking mutations data...\n")
  missing_cols <- setdiff(required_mut_cols, names(incommon_data$mutations))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in mutations: ", 
         paste(missing_cols, collapse = ", "))
  }
  cat("  ✓ All required columns present\n")
  
  # Check for NA values
  for (col in required_mut_cols) {
    na_count <- sum(is.na(incommon_data$mutations[[col]]))
    if (na_count > 0) {
      warning("  ⚠ ", na_count, " NA values in ", col)
    } else {
      cat("  ✓", col, "complete\n")
    }
  }
  
  # Check value ranges
  cat("\nChecking value ranges...\n")
  
  if (any(incommon_data$mutations$t_depth <= 0, na.rm = TRUE)) {
    warning("  ⚠ Some t_depth values are ≤ 0")
  } else {
    cat("  ✓ t_depth values > 0\n")
  }
  
  if (any(incommon_data$mutations$VAF < 0 | incommon_data$mutations$VAF > 1, na.rm = TRUE)) {
    warning("  ⚠ Some VAF values outside [0,1] range")
  } else {
    cat("  ✓ VAF values in valid range [0,1]\n")
  }
  
  if (any(incommon_data$mutations$PURITY < 0 | incommon_data$mutations$PURITY > 1, na.rm = TRUE)) {
    warning("  ⚠ Some PURITY values outside [0,1] range")
  } else {
    cat("  ✓ PURITY values in valid range [0,1]\n")
  }
  
  # Summary statistics
  cat("\nSummary Statistics:\n")
  cat("  Genes:", length(unique(incommon_data$mutations$Hugo_Symbol)), "\n")
  cat("  Samples:", incommon_data$parameters$n_samples, "\n")
  cat("  Patients:", incommon_data$parameters$n_patients, "\n")
  cat("  Mutations:", incommon_data$parameters$n_mutations, "\n")
  cat("  Mean VAF:", round(mean(incommon_data$mutations$VAF, na.rm = TRUE), 3), "\n")
  cat("  Mean t_depth:", round(mean(incommon_data$mutations$t_depth, na.rm = TRUE), 1), "\n")
  cat("  Mean purity:", round(mean(incommon_data$purity$PURITY), 3), "(fraction)\n")
  
  cat("\n✓ Data verification complete\n")
  invisible(TRUE)
}

#' Add copy number data to mutations for INCOMMON analysis
#'
#' @param incommon_data Object from prepare_incommon_data()
#' @return incommon_ready object with copy number added to mutations
#' @export
add_copy_number_to_mutations <- function(incommon_data) {
  
  if (!inherits(incommon_data, "incommon_ready")) {
    stop("Input must be an incommon_ready object from prepare_incommon_data()")
  }
  
  if (is.null(incommon_data$cna)) {
    stop("No CNA data available. Re-run prepare_incommon_data() to include CNA.")
  }
  
  mutations <- incommon_data$mutations
  cna_data <- incommon_data$cna
  
  cat("Matching copy numbers to mutations...\n")
  mutations$copy_number <- NA
  
  for (i in 1:nrow(mutations)) {
    gene <- mutations$Hugo_Symbol[i]
    sample <- mutations$SAMPLE_ID[i]
    
    if (gene %in% cna_data$GENE && sample %in% names(cna_data)) {
      mutations$copy_number[i] <- cna_data[cna_data$GENE == gene, sample]
    }
  }
  
  matched <- sum(!is.na(mutations$copy_number))
  cat("  Matched", matched, "of", nrow(mutations), "mutations\n")
  
  # Convert log2 ratio to absolute copy numbers
  # Formula: absolute_cn = 2 * 2^(log2_ratio)
  mutations$absolute_cn <- round(2 * 2^mutations$copy_number)
  
  # k_total for INCOMMON (ensure non-negative)
  mutations$k_total <- pmax(mutations$absolute_cn, 0, na.rm = TRUE)
  
  incommon_data$mutations <- mutations
  incommon_data$parameters$has_copy_number <- TRUE
  
  return(incommon_data)
}

#' Estimate eta parameter for INCOMMON
#'
#' Estimates the read rate per chromosome copy (eta) from diploid regions
#'
#' @param mutations Data frame with mutation data including k_total, VAF, t_depth, PURITY
#' @return Numeric value of estimated eta
#' @export
estimate_eta <- function(mutations) {
  
  cat("Estimating eta (reads per chromosome copy)...\n")
  
  if (!"k_total" %in% names(mutations)) {
    warning("No copy number data (k_total). Using fallback: eta = mean(depth)/2")
    eta_est <- mean(mutations$t_depth, na.rm = TRUE) / 2
    cat("  Fallback eta:", round(eta_est, 2), "\n")
    return(eta_est)
  }
  
  # Use diploid regions (k=2) with moderate VAF (likely heterozygous)
  diploid_muts <- mutations %>%
    dplyr::filter(!is.na(k_total), k_total == 2, 
                  VAF > 0.2, VAF < 0.6, 
                  !is.na(PURITY))
  
  if (nrow(diploid_muts) > 0) {
    cat("  Using", nrow(diploid_muts), "diploid mutations for estimation\n")
    # For diploid regions: expected depth ≈ 2*eta
    eta_est <- mean(diploid_muts$t_depth, na.rm = TRUE) / 2
    cat("  Estimated eta:", round(eta_est, 2), "\n")
    return(eta_est)
  }
  
  # Fallback
  warning("No suitable diploid mutations found. Using fallback method.")
  eta_est <- mean(mutations$t_depth, na.rm = TRUE) / 2
  cat("  Fallback eta:", round(eta_est, 2), "\n")
  return(eta_est)
}

#' Check INCOMMON data readiness
#'
#' Verifies that all required parameters for INCOMMON analysis are present
#'
#' @param incommon_data Object from prepare_incommon_data()
#' @param eta Optional: eta value (reads per chromosome). If NULL, will note it's needed.
#' @return Invisible list with readiness status
#' @export
check_incommon_readiness <- function(incommon_data, eta = NULL) {
  
  cat("\nCHECKING INCOMMON READINESS\n")
  cat("===========================\n\n")
  
  muts <- incommon_data$mutations
  
  # Check basic columns
  required <- c("Hugo_Symbol", "t_alt_count", "t_depth", "PURITY")
  missing <- setdiff(required, names(muts))
  
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  cat("✓ Basic columns present\n")
  
  # Check copy number
  has_cn <- "k_total" %in% names(muts)
  if (!has_cn) {
    cat("✗ Missing k_total (copy number)\n")
    cat("  Run: data <- add_copy_number_to_mutations(data)\n")
  } else {
    na_cn <- sum(is.na(muts$k_total))
    cat("✓ Copy number (k_total) available\n")
    cat("  Range:", paste(range(muts$k_total, na.rm = TRUE), collapse = " - "), "\n")
    if (na_cn > 0) {
      cat("  Warning:", na_cn, "mutations without copy number\n")
    }
  }
  
  # Check eta
  has_eta <- !is.null(eta)
  if (!has_eta) {
    cat("✗ eta (reads per chromosome) not provided\n")
    cat("  Run: eta <- estimate_eta(data$mutations)\n")
  } else {
    cat("✓ eta =", round(eta, 2), "\n")
  }
  
  # Parameter summary
  cat("\n=== INCOMMON Parameters ===\n")
  cat("For each mutation i, you have:\n")
  cat("  r_i (variant reads):", paste(range(muts$t_alt_count), collapse = " "), "\n")
  cat("  d_i (total depth):  ", paste(range(muts$t_depth), collapse = " "), "\n")
  cat("  π (purity):         ", paste(round(range(muts$PURITY), 3), collapse = " - "), 
      "(fraction)\n")
  
  if (has_cn) {
    cat("  k_i (copy number):  ", paste(range(muts$k_total, na.rm = TRUE), collapse = " "), "\n")
  }
  
  if (has_eta) {
    cat("  η (reads per chrom):", round(eta, 2), "\n")
  }
  
  cat("\nDataset summary:\n")
  cat("  Mutations:", nrow(muts), "\n")
  cat("  Samples:", length(unique(muts$SAMPLE_ID)), "\n")
  cat("  Patients:", length(unique(muts$PATIENT_ID)), "\n")
  cat("  Genes:", length(unique(muts$Hugo_Symbol)), "\n")
  
  ready <- has_cn && has_eta
  
  cat("\n")
  if (ready) {
    cat("✓ READY FOR INCOMMON ANALYSIS\n")
  } else {
    cat("✗ NOT READY - complete steps above\n")
  }
  
  return(invisible(list(ready = ready, eta = eta, has_cn = has_cn)))
}
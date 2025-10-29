#' @importFrom magrittr %>%
NULL

#' Download MSK-CHORD data directly from cBioPortal API
#'
#' This function completely bypasses cBioPortalData and downloads the data
#' directly from the cBioPortal REST API, then building a MultiAssayExperiment
#'
#' @param study Study identifier (default: "msk_chord_2024")
#' @param api_url Base URL for cBioPortal API
#' @param validate_data Whether to run data validation checks
#' @return MultiAssayExperiment object containing MSK-CHORD data
#' @export
download_msk_chord <- function(study = "msk_chord_2024",
                               api_url = "https://www.cbioportal.org/api",
                               validate_data = TRUE,
                               cache_dir = NULL,
                               force_refresh = FALSE) {
  
  if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
    stop("Package 'cBioPortalData' required. Install with: BiocManager::install('cBioPortalData')", 
         call. = FALSE)
  }
  
  if (!is.character(study) || length(study) != 1) {
    stop("'study' must be a single character string", call. = FALSE)
  }
  
  # Define cache directory
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("MSK_CHORD", "cache")
  }
  
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  tarball_path <- file.path(cache_dir, paste0(study, ".tar.gz"))
  extract_dir <- file.path(cache_dir, "extracted")
  study_dir <- file.path(extract_dir, study)
  
  # Check if already exists in cache
  if (!force_refresh && dir.exists(study_dir)) {
    message("Loading MSK-CHORD data from cache...")
    
    tryCatch({
      data <- cBioPortalData::loadStudy(filepath = study_dir)
      
      attr(data, "study_id") <- study
      attr(data, "download_date") <- file.info(study_dir)$mtime
      attr(data, "source") <- "cache"
      
      if (validate_data) {
        validate_msk_chord_data(data)
      }
      
      message("Successfully loaded MSK-CHORD data from cache")
      message(sprintf("  - Cache location: %s", study_dir))
      return(data)
      
    }, error = function(e) {
      message("Cache corrupted or invalid, re-downloading...")
      message("  Error: ", e$message)
      # Continue with download
    })
  }
  
  # Download and load dataset
  tryCatch({
    message("Downloading MSK-CHORD tarball from cBioPortal...")
    
    # URL of dataset hub
    url <- sprintf(
      "https://cbioportal-datahub.s3.amazonaws.com/%s.tar.gz",
      study
    )
    
    # Download with progress
    message(sprintf("  - Downloading from: %s", url))
    download.file(url, destfile = tarball_path, mode = "wb", quiet = FALSE)
    message("  - Downloaded tarball successfully")
    
    # Extract using untar
    message("  - Extracting files...")
    utils::untar(
      tarfile = tarball_path,
      exdir = extract_dir
    )
    message("  - Files extracted to:", extract_dir)
    
    # Load as MultiAssayExperiment using loadStudy
    message("  - Building MultiAssayExperiment...")
    data <- cBioPortalData::loadStudy(filepath = study_dir)
    
    # Remove the tarball to save space
    if (file.exists(tarball_path)) {
      file.remove(tarball_path)
      message("  - Removed tarball to save space")
    }
    
    attr(data, "study_id") <- study
    attr(data, "download_date") <- Sys.Date()
    attr(data, "source") <- "download"
    
    if (validate_data) {
      validate_msk_chord_data(data)
    }
    
    message("Successfully loaded MSK-CHORD data")
    message(sprintf("  - Samples: %d", ncol(data)))
    message(sprintf("  - Experiments: %d", length(data)))
    message(sprintf("  - Assays: %s", paste(names(data), collapse = ", ")))
    message(sprintf("  - Cached in: %s", study_dir))
    
    return(data)
    
  }, error = function(e) {
    stop("Failed to download MSK-CHORD data: ", e$message,
         "\n\nTroubleshooting:",
         "\n1. Check if URL is correct: https://cbioportal-datahub.s3.amazonaws.com/",
         "\n2. Try downloading manually from: https://www.cbioportal.org/datasets",
         "\n3. Check available disk space",
         "\n4. Verify internet connection",
         call. = FALSE)
  })
}

#' Validate MSK-CHORD data structure
#' 
#' @param data MultiAssayExperiment object
#' @return Invisible validation results
#' @keywords internal
validate_msk_chord_data <- function(data) {
  
  if (!inherits(data, "MultiAssayExperiment")) {
    stop("Data is not a MultiAssayExperiment object", call. = FALSE)
  }
  
  # Robust check for number of columns (samples)
  num_cols <- try(ncol(data), silent = TRUE)
  if (inherits(num_cols, "try-error") || length(num_cols) == 0 || num_cols == 0) {
    warning("MultiAssayExperiment has no samples", call. = FALSE)
  }
  
  # Robust check for number of experiments
  num_exps <- try(length(data), silent = TRUE)
  if (inherits(num_exps, "try-error") || length(num_exps) == 0 || num_exps == 0) {
    warning("MultiAssayExperiment has no experiments", call. = FALSE)
  }
  
  message("Data validation complete")
  invisible(TRUE)
}
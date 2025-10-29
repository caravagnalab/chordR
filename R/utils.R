#' @importFrom magrittr %>%
NULL

#' Pipe operator
#' 
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder
#' @param rhs A function call using the magrittr semantics
#' @return The result of calling rhs(lhs)
NULL

#' Null coalescing operator
#' 
#' @param x First value
#' @param y Second value (default if x is NULL)
#' @return x if not NULL, otherwise y
#' @keywords internal
#' @name null-coalesce
#' @rdname null-coalesce
NULL

#' @rdname null-coalesce
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

#' Get unique patient IDs from timeline object
#' 
#' @param timeline_obj MSKTimeline object
#' @return Character vector of unique patient IDs
#' @export
#' @examples
#' \dontrun{
#' # Create example timeline
#' tl <- new("MSKTimeline", 
#'           data = list(labs = data.frame(PATIENT_ID = c("P1", "P2"))),
#'           n_patients = 2L)
#' get_patient_ids(tl)
#' }
get_patient_ids <- function(timeline_obj) {
  
  if (!inherits(timeline_obj, "MSKTimeline")) {
    stop("Input must be an MSKTimeline object", call. = FALSE)
  }
  
  if (length(timeline_obj@data) == 0) {
    return(character(0))
  }
  
  unique(unlist(lapply(timeline_obj@data, function(x) {
    if ("PATIENT_ID" %in% names(x)) {
      x$PATIENT_ID
    } else {
      character(0)
    }
  })))
}

#' Check if required packages are available
#' 
#' @param packages Character vector of package names
#' @param purpose Description of what packages are needed for
#' @return Invisible logical vector
#' @keywords internal
check_required_packages <- function(packages, purpose = "this functionality") {
  
  available <- vapply(packages, requireNamespace, 
                      logical(1), quietly = TRUE)
  
  if (!all(available)) {
    missing <- packages[!available]
    stop("Required package(s) for ", purpose, ": ", 
         paste(missing, collapse = ", "), "\n",
         "Install with: install.packages(c('", 
         paste(missing, collapse = "', '"), "'))",
         call. = FALSE)
  }
  
  invisible(available)
}

#' Combine data frames with different columns
#' 
#' @param df_list List of data frames to combine
#' @return Combined data frame with all columns
#' @export
#' @examples
#' \dontrun{
#' df1 <- data.frame(a = 1:3, b = 4:6)
#' df2 <- data.frame(a = 7:9, c = 10:12)
#' bind_rows_safe(list(df1, df2))
#' }
bind_rows_safe <- function(df_list) {
  
  if (length(df_list) == 0) {
    return(data.frame())
  }
  
  if (length(df_list) == 1) {
    return(df_list[[1]])
  }
  
  # Get all column names
  all_colnames <- lapply(df_list, names)
  all_cols <- unique(unlist(all_colnames))
  
  # Align all data frames
  aligned_list <- lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    
    # Add missing columns with NA
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    
    # Reorder to match
    df[, all_cols, drop = FALSE]
  })
  
  # Combine
  result <- do.call(rbind, aligned_list)
  rownames(result) <- NULL
  
  return(result)
}

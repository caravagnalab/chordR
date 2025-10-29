#' @importFrom methods new setClass setMethod show
NULL

#' MSK-CHORD Timeline Data Class
#' 
#' An S4 class to represent timeline data from MSK-CHORD dataset
#' 
#' @slot data Named list of timeline data frames
#' @slot metadata List of metadata information
#' @slot cancer_types Character vector of included cancer types
#' @slot n_patients Integer number of unique patients
#' 
#' @exportClass MSKTimeline
setClass("MSKTimeline",
         slots = list(
           data = "list",
           metadata = "list", 
           cancer_types = "character",
           n_patients = "integer"
         ),
         prototype = list(
           data = list(),
           metadata = list(),
           cancer_types = character(0),
           n_patients = 0L
         )
)

#' Show method for MSKTimeline
#' 
#' @param object An MSKTimeline object
#' @export
setMethod("show", "MSKTimeline", function(object) {
  cat("MSK-CHORD Timeline Object\n")
  cat("========================\n")
  cat("Patients:", object@n_patients, "\n")
  cat("Cancer Types:", paste(object@cancer_types, collapse = ", "), "\n")
  cat("Timeline Types:", length(object@data), "\n")
  
  if (length(object@data) > 0) {
    total_events <- sum(sapply(object@data, nrow))
    cat("Total Events:", total_events, "\n")
  }
})

#' Summary method for MSKTimeline
#' 
#' @param object An MSKTimeline object
#' @export
setMethod("summary", "MSKTimeline", function(object) {
  list(
    n_patients = object@n_patients,
    n_events = sum(sapply(object@data, nrow)),
    n_timeline_types = length(object@data),
    cancer_types = object@cancer_types
  )
})
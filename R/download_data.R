#' Download data from cBioPortal
#'
#' @param study The identifier of the study on cBioPortal
#' @return A MultiAssayExperiment object
#' @export
#' @importFrom cBioPortalData cBioDataPack

download_data = function(study = "msk_chord_2024"){

  msk_chord = cBioPortalData::cBioDataPack(study, ask = FALSE)
}

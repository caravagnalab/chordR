#' Subset
#'
#' @param x The input MultiAssayExperiment object
#' @return A MultiAssayExperiment object
#' @export
#' @importFrom cBioPortalData cBioDataPack

subset_data = function(x){

  msk_chord = cBioPortalData::cBioDataPack(study, ask = FALSE)
}

# xx = subsetByColData(x = msk_chord, y = "P-0002869")
# xx@ExperimentList@listData$mutations@assays %>% names()
# 
# msk_chord@ExperimentList@listData$mutations@assays %>% class()
# 
# yy = subsetByColumn(x = msk_chord, y = "P-0002869")
# yy@ExperimentList@listData$mutations@assays %>% as_tibble()
# 
# zz = subsetByRow(x = msk_chord, y = "P-0002869")
# zz@ExperimentList@listData$mutations@assays %>% as_tibble()
# 
# colData(x = msk_chord) %>% colnames()

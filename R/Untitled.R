experiments(x = msk_chord)


sampleMap(x = msk_chord) %>% as_tibble()

cn = colnames(colData(x = msk_chord))

# Mutations
msk_chord@ExperimentList@listData$mutations@assays %>% as_tibble()

# Clinical data
colData(x = msk_chord) %>% as_tibble()

# Structural Variants
msk_chord@metadata$sv %>% as_tibble()

# Timelines
msk_chord@metadata$`timeline_ca_15-3_labs` %>% as_tibble()

msk_chord@metadata$timeline_psa_labs %>%
  as_tibble() %>%
  dplyr::slice_head(n = 500) %>%
  ggplot2::ggplot(ggplot2::aes(x = START_DATE, y = log10(RESULT)))+
  ggplot2::geom_line(ggplot2::aes(group = PATIENT_ID, color = PATIENT_ID))+
  ggplot2::geom_point(ggplot2::aes(group = PATIENT_ID, color = PATIENT_ID))

msk_chord@metadata$timeline_gleason %>%
  as_tibble() %>%
  dplyr::slice_head(n = 100) %>%
  dplyr::mutate(gleason = ifelse(GLEASON_SCORE >=7, 'high', 'low')) %>%
  ggplot2::ggplot(ggplot2::aes(x = START_DATE, y = (GLEASON_SCORE)))+
  ggplot2::geom_line(ggplot2::aes(group = PATIENT_ID, color = PATIENT_ID))+
  ggplot2::geom_point(ggplot2::aes(group = PATIENT_ID, color = PATIENT_ID))+
  ggplot2::ylim(0,100)+ggplot2::xlim(0,1000)+
  ggplot2::facet_wrap(~gleason)

xx@metadata$timeline_progression %>% as_tibble()

msk_chord@metadata %>% as_tibble()

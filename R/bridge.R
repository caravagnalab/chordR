#' Convert chordR output to INCOMMON input format
#'
#' @param incommon_data Output from chordR::prepare_incommon_data() + add_copy_number_to_mutations()
#' @param cohort Output from chordR::extract_by_cancer_type()
#' @return List with genomic_data and clinical_data ready for INCOMMON::init()
#' @export
chordr_to_incommon <- function(incommon_data, cohort) {
  
  cat("=== CONVERSION chordR → INCOMMON ===\n\n")
  
  mutations <- incommon_data$mutations
  
  # ==========================================
  # 1. PREPARE GENOMIC_DATA
  # ==========================================
  
  cat("Step 1: Prepare genomic_data...\n")
  
  genomic_data <- mutations %>%
    dplyr::transmute(
      # === Mandatory Columns ===
      sample = SAMPLE_ID,
      gene = Hugo_Symbol,
      NV = as.integer(t_alt_count),
      DP = as.integer(t_depth),
      VAF = VAF,
      
      # === Genomic Position ===
      chr = as.character(seqnames),
      from = as.integer(start),
      to = as.integer(end),
      
      # === Alleles ===
      ref = if("Reference_Allele" %in% names(mutations)) {
        as.character(Reference_Allele)
      } else {
        NA_character_
      },
      
      alt = if("Tumor_Seq_Allele2" %in% names(mutations)) {
        as.character(Tumor_Seq_Allele2)
      } else {
        NA_character_
      },
      
      # === Protein Annotation ===
      HGVSp_Short = if("HGVSp_Short" %in% names(mutations)) {
        HGVSp_Short
      } else {
        NA_character_
      },
      
      # === Extra info ===
      variant_classification = if("Variant_Classification" %in% names(mutations)) {
        Variant_Classification
      } else {
        NA_character_
      },
      
      variant_type = if("Variant_Type" %in% names(mutations)) {
        Variant_Type
      } else {
        NA_character_
      },
      
      consequence = if("Consequence" %in% names(mutations)) {
        Consequence
      } else {
        NA_character_
      }
    )
  
  # Quality filtering
  initial_rows <- nrow(genomic_data)
  
  genomic_data <- genomic_data %>%
    dplyr::filter(
      !is.na(sample),
      !is.na(gene),
      !is.na(NV),
      !is.na(DP),
      !is.na(VAF),
      DP > 0,
      NV >= 0,
      NV <= DP,
      VAF >= 0,
      VAF <= 1
    )
  
  removed <- initial_rows - nrow(genomic_data)
  if (removed > 0) {
    cat("  ⚠️  Removed", removed, "record with null/invalid data\n")
  }
  
  cat("  ✓ genomic_data:", nrow(genomic_data), "mutations\n")
  cat("  ✓ Samples:", length(unique(genomic_data$sample)), "\n")
  cat("  ✓ Genes:", length(unique(genomic_data$gene)), "\n\n")
  
  # ==========================================
  # 2. CLINICAL_DATA
  # ==========================================
  
  cat("Step 2: Prepare clinical_data...\n")
  
  # Extract purity from mutations
  purity_from_muts <- incommon_data$purity %>%
    dplyr::transmute(
      SAMPLE_ID = SAMPLE_ID,
      purity_from_mutations = PURITY
    ) %>%
    dplyr::distinct()
  
  # Merge with complete clinical data
  if (!is.null(cohort$clinical)) {
    
    clinical_data <- cohort$clinical %>%
      dplyr::transmute(
        # === Identifiers ===
        sample = SAMPLE_ID,
        patient = PATIENT_ID,
        
        # === PURITY  ===
        purity = if("TUMOR_PURITY" %in% names(.)) {
          TUMOR_PURITY / 100 # fraction 
        } else {
          NA_real_
        },
        
        # === TUMOR TYPE  ===
        tumor_type = if("ONCOTREE_CODE" %in% names(.)) {
          ONCOTREE_CODE
        } else if("CANCER_TYPE" %in% names(.)) {
          CANCER_TYPE
        } else {
          cohort$cancer_type
        },
        
        cancer_type = if("CANCER_TYPE" %in% names(.)) {
          CANCER_TYPE
        } else {
          NA_character_
        },
        
        cancer_type_detailed = if("CANCER_TYPE_DETAILED" %in% names(.)) {
          CANCER_TYPE_DETAILED
        } else {
          NA_character_
        },
        
        # === SAMPLE CHARACTERISTICS ===
        SAMPLE_TYPE = if("SAMPLE_TYPE" %in% names(.)) {
          SAMPLE_TYPE
        } else {
          NA_character_
        },
        
        SAMPLE_CLASS = if("SAMPLE_CLASS" %in% names(.)) {
          SAMPLE_CLASS
        } else {
          NA_character_
        },
        
        PRIMARY_SITE = if("PRIMARY_SITE" %in% names(.)) {
          PRIMARY_SITE
        } else {
          NA_character_
        },
        
        METASTATIC_SITE = if("METASTATIC_SITE" %in% names(.)) {
          METASTATIC_SITE
        } else {
          NA_character_
        },
        
        # === METASTATIC SITES ===
        met_adrenal = if("ADRENAL_GLANDS" %in% names(.)) {
          ADRENAL_GLANDS
        } else {
          NA_integer_
        },
        
        met_bone = if("BONE" %in% names(.)) {
          BONE
        } else {
          NA_integer_
        },
        
        met_brain = if("CNS_BRAIN" %in% names(.)) {
          CNS_BRAIN
        } else {
          NA_integer_
        },
        
        met_liver = if("LIVER" %in% names(.)) {
          LIVER
        } else {
          NA_integer_
        },
        
        met_lung = if("LUNG" %in% names(.)) {
          LUNG
        } else {
          NA_integer_
        },
        
        met_lymph_nodes = if("LYMPH_NODES" %in% names(.)) {
          LYMPH_NODES
        } else {
          NA_integer_
        },
        
        # MET_COUNT
        MET_COUNT = rowSums(dplyr::select(., dplyr::starts_with("met_")), na.rm = TRUE),
        
        # === SURVIVAL DATA ===
        OS_STATUS = if("OS_STATUS" %in% names(.)) {
          OS_STATUS
        } else {
          NA_character_
        },
        
        OS_MONTHS = if("OS_MONTHS" %in% names(.)) {
          OS_MONTHS
        } else {
          NA_real_
        },
        
        # === DEMOGRAPHICS ===
        age = if("CURRENT_AGE_DEID" %in% names(.)) {
          CURRENT_AGE_DEID
        } else {
          NA_real_
        },
        
        sex = if("GENDER" %in% names(.)) {
          GENDER
        } else {
          NA_character_
        },
        
        race = if("RACE" %in% names(.)) {
          RACE
        } else {
          NA_character_
        },
        
        ethnicity = if("ETHNICITY" %in% names(.)) {
          ETHNICITY
        } else {
          NA_character_
        },
        
        # === STAGING ===
        stage = if("STAGE_HIGHEST_RECORDED" %in% names(.)) {
          STAGE_HIGHEST_RECORDED
        } else {
          NA_character_
        },
        
        # === MOLECULAR FEATURES ===
        MSI_TYPE = if("MSI_TYPE" %in% names(.)) {
          MSI_TYPE
        } else {
          NA_character_
        },
        
        MSI_SCORE = if("MSI_SCORE" %in% names(.)) {
          MSI_SCORE
        } else {
          NA_real_
        },
        
        TMB_NONSYNONYMOUS = if("TMB_NONSYNONYMOUS" %in% names(.)) {
          TMB_NONSYNONYMOUS
        } else {
          NA_real_
        },
        
        PDL1_POSITIVE = if("PDL1_POSITIVE" %in% names(.)) {
          PDL1_POSITIVE
        } else {
          NA_character_
        },
        
        # === BIOMARKERS ===
        HR = if("HR" %in% names(.)) {
          HR
        } else {
          NA_character_
        },
        
        HER2 = if("HER2" %in% names(.)) {
          HER2
        } else {
          NA_character_
        },
        
        # === LUNG CANCER SPECIFIC ===
        smoking_status = if("SMOKING_PREDICTIONS_3_CLASSES" %in% names(.)) {
          SMOKING_PREDICTIONS_3_CLASSES
        } else {
          NA_character_
        },
        
        # === SEQUENCING INFO ===
        gene_panel = if("GENE_PANEL" %in% names(.)) {
          GENE_PANEL
        } else {
          NA_character_
        },
        
        sample_coverage = if("SAMPLE_COVERAGE" %in% names(.)) {
          SAMPLE_COVERAGE
        } else {
          NA_real_
        }
      ) %>%
      dplyr::distinct()
    
    # Merge purity from mutations if missing
    clinical_data <- clinical_data %>%
      dplyr::left_join(purity_from_muts, by = c("sample" = "SAMPLE_ID")) %>%
      dplyr::mutate(
        purity = dplyr::coalesce(purity, purity_from_mutations)
      ) %>%
      dplyr::select(-purity_from_mutations)
    
  } else {
    # Only purity if clinical not available
    clinical_data <- purity_from_muts %>%
      dplyr::rename(sample = SAMPLE_ID, purity = purity_from_mutations)
    
    if (!is.null(cohort$cancer_type)) {
      clinical_data$tumor_type <- cohort$cancer_type
    }
  }

  # verify purity is in fraction format
  if (any(clinical_data$purity > 1, na.rn = TRUE)) {
    warning("⚠️ - Some purity values > 1 detected. Converting from percentage to fraction")
    clinical_data$purity <- clinical_data$purity / 100
  }
  
  # TMB from mutations
  tmb_calculated <- genomic_data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(tmb = n(), .groups = "drop")
  
  clinical_data <- clinical_data %>%
    dplyr::left_join(tmb_calculated, by = "sample")
  
  # Remove column patient
  if ("patient" %in% names(clinical_data)) {
    clinical_data <- clinical_data %>% dplyr::select(-patient)
  }
  
  cat("  ✓ clinical_data:", nrow(clinical_data), "samples\n")
  cat("  ✓ Purity range:", 
      paste(round(range(clinical_data$purity, na.rm = TRUE), 3), 
            collapse = " - "), "(fraction)\n") #### CHANGED!!!
  
  # Report avilable fields
  available_fields <- sum(!sapply(clinical_data[1,], is.na))
  cat("  ✓ Available fields:", available_fields, "di", ncol(clinical_data), "\n")
  
  # verify key data availability
  if ("OS_STATUS" %in% names(clinical_data)) {
    n_survival <- sum(!is.na(clinical_data$OS_STATUS))
    cat("  ℹ️  Survival data:", n_survival, "samples\n")
  }
  
  if ("SAMPLE_TYPE" %in% names(clinical_data)) {
    sample_types <- table(clinical_data$SAMPLE_TYPE, useNA = "ifany")
    cat("  ℹ️  Sample types:", paste(names(sample_types), sample_types, 
                                     sep = "=", collapse = ", "), "\n")
  }
  
  if ("MET_COUNT" %in% names(clinical_data)) {
    n_with_mets <- sum(clinical_data$MET_COUNT > 0, na.rm = TRUE)
    cat("  ℹ️  Samples with metastases:", n_with_mets, "\n")
  }
  
  cat("\n")
  
  # ==========================================
  # 3. Verify Compatibility
  # ==========================================
  
  cat("Step 3: Verify compatibility...\n")
  
  samples_genomic <- unique(genomic_data$sample)
  samples_clinical <- unique(clinical_data$sample)
  
  only_genomic <- setdiff(samples_genomic, samples_clinical)
  only_clinical <- setdiff(samples_clinical, samples_genomic)
  matching <- intersect(samples_genomic, samples_clinical)
  
  if (length(only_genomic) > 0) {
    cat("  ⚠️ ", length(only_genomic), 
        "samples in genomic_data without clinical_data\n")
  }
  
  if (length(only_clinical) > 0) {
    cat("  ⚠️ ", length(only_clinical), 
        "samples in clinical_data without genomic_data\n")
  }
  
  cat("  ✓ Samples matching:", length(matching), "\n\n")
  
  # ==========================================
  # 4. QUALITY CHECKS
  # ==========================================
  
  cat("Step 4: Quality checks...\n")
  
  checks <- list()
  
  # Check VAF
  invalid_vaf <- sum(genomic_data$VAF < 0 | genomic_data$VAF > 1, na.rm = TRUE)
  checks$vaf <- invalid_vaf == 0
  if (invalid_vaf == 0) {
    cat("  ✓ VAF in valid range\n")
  } else {
    cat("  ⚠️ ", invalid_vaf, "mutations with VAF out of range [0,1]\n")
  }
  
  # Check read counts
  invalid_counts <- sum(genomic_data$NV > genomic_data$DP, na.rm = TRUE)
  checks$counts <- invalid_counts == 0
  if (invalid_counts == 0) {
    cat("  ✓ Valid Read Counts (NV ≤ DP)\n")
  } else {
    cat("  ⚠️ ", invalid_counts, "mutations with NV > DP\n")
  }
  
  # Check purity
  invalid_purity <- sum(clinical_data$purity < 0 | 
                          clinical_data$purity > 1, na.rm = TRUE)
  checks$purity <- invalid_purity == 0
  if (invalid_purity == 0) {
    cat("  ✓ Purity in valid range [0,1]\n")
  } else {
    cat("  ⚠️ ", invalid_purity, "samples with purity out of range [0,1]\n")
  }
  
  # Check genomic position
  missing_position <- sum(is.na(genomic_data$chr) | is.na(genomic_data$from))
  checks$position <- missing_position == 0
  if (missing_position == 0) {
    cat("  ✓ Complete Genomic Position\n")
  } else {
    cat("  ⚠️ ", missing_position, 
        "mutations without genomic positions (chr/from)\n")
  }
  
  # Check alleles 
  missing_alleles <- sum(is.na(genomic_data$ref) | is.na(genomic_data$alt))
  if (missing_alleles == 0) {
    cat("  ✓ Complete Alleles (ref/alt)\n")
  } else {
    cat("  ℹ️ ", missing_alleles, 
        "mutations without alleles (ref/alt)\n")
  }
  
  # Check tumor_type
  if ("tumor_type" %in% names(clinical_data)) {
    missing_tumor_type <- sum(is.na(clinical_data$tumor_type))
    if (missing_tumor_type == 0) {
      cat("  ✓ tumor_type available for all samples\n")
      tumor_types <- unique(clinical_data$tumor_type)
      if (length(tumor_types) <= 3) {
        cat("    Tumor types:", paste(tumor_types, collapse = ", "), "\n")
      }
    } else {
      cat("  ℹ️ ", missing_tumor_type, "samples without tumor_type\n")
      cat("      INCOMMON will use priors PANCA (pan-cancer)\n")
    }
  } else {
    cat("  ℹ️  tumor_type not available\n")
    cat("      INCOMMON will use priors PANCA (pan-cancer)\n")
  }
  
  cat("\n")
  
  # ==========================================
  # 5. SUMMARY
  # ==========================================
  
  cat("=== SUMMARY ===\n")
  cat("genomic_data:\n")
  cat("  - Mutations:", nrow(genomic_data), "\n")
  cat("  - Samples:", length(unique(genomic_data$sample)), "\n")
  cat("  - Genes:", length(unique(genomic_data$gene)), "\n")
  cat("  - Mean DP:", round(mean(genomic_data$DP), 1), "\n")
  cat("  - Mean VAF:", round(mean(genomic_data$VAF), 3), "\n")
  
  cat("\nclinical_data:\n")
  cat("  - Samples:", nrow(clinical_data), "\n")
  cat("  - Mean purity:", round(mean(clinical_data$purity, na.rm = TRUE), 3), "\n")
  cat("  - Purity range:", paste(round(range(clinical_data$purity, na.rm=TRUE), 3), collapse = " - "), "\n")
  cat("  - Median TMB:", median(clinical_data$tmb, na.rm = TRUE), "\n")
  
  # Analysis for tumor type
  if ("tumor_type" %in% names(clinical_data)) {
    tumor_summary <- clinical_data %>%
      dplyr::group_by(tumor_type) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(n))
    
    if (nrow(tumor_summary) <= 5) {
      cat("  - Tumor types:\n")
      for (i in 1:nrow(tumor_summary)) {
        cat("    •", tumor_summary$tumor_type[i], ":", 
            tumor_summary$n[i], "samples\n")
      }
    }
  }
  
  all_checks_pass <- all(unlist(checks))
  
  cat("\n")
  if (all_checks_pass) {
    cat("✓ DATA READY FOR INCOMMON!\n")
  } else {
    cat("⚠️  Some checks failed - review warnings above\n")
  }
  
  cat("\nNext steps:\n")
  cat("1. library(INCOMMON)\n")
  cat("2. data('priors_pcawg_hmf', package = 'INCOMMON')\n")
  cat("3. data('priors_eta', package = 'INCOMMON')\n")
  cat("4. x <- INCOMMON::init(genomic_data, clinical_data)\n")
  cat("5. x <- INCOMMON::classify(x, priors_k_m, priors_eta)\n")
  
  # ==========================================
  # RETURN
  # ==========================================
  
  result <- list(
    genomic_data = genomic_data,
    clinical_data = clinical_data,
    conversion_summary = list(
      n_mutations = nrow(genomic_data),
      n_samples_genomic = length(unique(genomic_data$sample)),
      n_samples_clinical = nrow(clinical_data),
      n_matching_samples = length(matching),
      checks_passed = checks,
      warnings = list(
        samples_only_genomic = only_genomic,
        samples_only_clinical = only_clinical,
        invalid_vaf = invalid_vaf,
        invalid_counts = invalid_counts,
        invalid_purity = invalid_purity,
        missing_position = missing_position,
        missing_alleles = missing_alleles
      ),
      available_analyses = list(
        survival = "OS_STATUS" %in% names(clinical_data) && 
          sum(!is.na(clinical_data$OS_STATUS)) > 0,
        metastatic_propensity = "SAMPLE_TYPE" %in% names(clinical_data) &&
          "MET_COUNT" %in% names(clinical_data),
        metastatic_tropism = "METASTATIC_SITE" %in% names(clinical_data) ||
          any(grepl("^met_", names(clinical_data)))
      )
    )
  )
  
  class(result) <- c("chordr_incommon_bridge", "list")
  
  return(result)
}


#' Print method for chordr_incommon_bridge
#' @export
print.chordr_incommon_bridge <- function(x, ...) {
  cat("chordR → INCOMMON Bridge Object\n")
  cat("===============================\n\n")
  
  cat("genomic_data:", nrow(x$genomic_data), "mutations,",
      length(unique(x$genomic_data$sample)), "samples\n")
  
  cat("clinical_data:", nrow(x$clinical_data), "samples\n")
  
  cat("\nQuality Checks:\n")
  checks <- x$conversion_summary$checks_passed
  for (check_name in names(checks)) {
    symbol <- if(checks[[check_name]]) "✓" else "✗"
    cat("  ", symbol, check_name, "\n")
  }
  
  cat("\nReady for: INCOMMON::init(genomic_data, clinical_data)\n")
  
  invisible(x)
}


#' Verify chordR data compatibility with INCOMMON
#'
#' @param incommon_data Output from prepare_incommon_data()
#' @return Logical indicating if data is compatible
#' @export
verify_incommon_compatibility <- function(incommon_data) {
  
  cat("VERIFY COMPATIBILITY chordR → INCOMMON\n")
  cat("========================================\n\n")
  
  mutations <- incommon_data$mutations
  
  # Required fields mapping
  required_mapping <- list(
    sample = c("SAMPLE_ID", "sample"),
    gene = c("Hugo_Symbol", "gene"),
    NV = c("t_alt_count", "NV"),
    DP = c("t_depth", "DP"),
    VAF = c("VAF")
  )
  
  all_ok <- TRUE
  
  for (incommon_col in names(required_mapping)) {
    possible_cols <- required_mapping[[incommon_col]]
    found <- any(possible_cols %in% names(mutations))
    
    if (found) {
      which_col <- possible_cols[possible_cols %in% names(mutations)][1]
      cat("✓", incommon_col, "←", which_col, "\n")
    } else {
      cat("✗", incommon_col, "← MANCANTE (opzioni:", 
          paste(possible_cols, collapse = ", "), ")\n")
      all_ok <- FALSE
    }
  }
  
  cat("\n")
  
  # Optional but useful fields
  optional_mapping <- list(
    chr = c("seqnames", "chr"),
    from = c("start", "from"),
    ref = c("Reference_Allele", "ref"),
    alt = c("Tumor_Seq_Allele2", "alt")
  )
  
  cat("Optional fields:\n")
  for (incommon_col in names(optional_mapping)) {
    possible_cols <- optional_mapping[[incommon_col]]
    found <- any(possible_cols %in% names(mutations))
    
    if (found) {
      which_col <- possible_cols[possible_cols %in% names(mutations)][1]
      cat("✓", incommon_col, "←", which_col, "\n")
    } else {
      cat("○", incommon_col, "← not found (options:", 
          paste(possible_cols, collapse = ", "), ")\n")
    }
  }
  
  cat("\n")
  
  # Check purity
  if (!is.null(incommon_data$purity)) {
    if ("PURITY" %in% names(incommon_data$purity)) {
      cat("✓ purity ← PURITY (clinical_data)\n")
    } else {
      cat("✗ PURITY not in incommon_data$purity\n")
      all_ok <- FALSE
    }
  } else {
    cat("✗ incommon_data$purity is NULL\n")
    all_ok <- FALSE
  }
  
  cat("\n")
  
  if (all_ok) {
    cat("✓ Compatible\n")
    cat("Use: chordr_to_incommon(incommon_data, cohort)\n")
  } else {
    cat("✗ Not compatible\n")
    cat("Fix problems above before continue.\n")
  }
  
  return(invisible(all_ok))
}
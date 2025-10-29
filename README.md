# chordR
# An R package for MSK-CHORD dataset

## Overview
The `chordR` package is designed to facilitate access to and analysis of the MSK-CHORD (Clinicogenomic, Harmonized Oncologic Real-World Dataset) from Memorial Sloan Kettering Cancer Center. This comprehensive dataset integrates clinical, genomic, and real-world data from 24,950 patients across five major cancer types, enabling advanced cancer genomics research and biomarker discovery.

## About MSK-CHORD
MSK-CHORD is a groundbreaking dataset that combines:
* **Natural Language Processing (NLP)** annotations from clinical notes and radiology reports
* **Tumor genomic profiling** via MSK-IMPACT targeted sequencing
* **Clinical data** including demographic, treatments, and outcomes
* **Timeline data** tracking patient events over time

### Dataset Composition
The dataset includes 24,950 patients across five principal cancer types:

| Cancer Type | Number of Patients |
|-------------|--------------------|
| Non-Small Cell Lung Cancer (NSCLC) | 7,809 |
| Breast Cancer | 5,368 |
| Colorectal Cancer | 5,543 |
| Prostate Cancer | 3,211 |
| Pancreatic Cancer | 3,109 |

### Key Features
* 705,241 radiology reports annotated with NLP
* Genomic alterations from MSK-IMPACT sequencing (341-gene panel)
* Survival outcomes and clinical follow-up
* Metastatic site information extracted from imaging reports

## Installation

```R
# install from GitHub
devtools::install_github("caravagnalab/chordR")

# load the package
library(chordR)
```

## Core Functions

### 1. Data Download and Access

`download_msk_chord()`

Download the complete MSK-CHORD dataset from cBioPortal.

```R
# basic download
msk_data <- download_msk_chord()

# with validation
msk_data <- download_msk_chord(
    study = "msk_chord_2024", # set as default
    validate_data = TRUE
)
```

**Parameters:**
* `study`: Study identifier (default: "msk_chord_2024")
* `validate_data`: Whether to run data validation checks (default: TRUE)
* `cache_dir`: Optional directory for caching downloaded data
* `force_refresh`: Force re-download even if cached data exists (default: FALSE)

**Returns:** A `MultAssayExperiment` object containing all MSK-CHORD data.

### 2. Timeline Data Extraction
`extract_msk_timeline()`

Extract and structure timeline data including lab results, treatments, and clinical events.

```R
# extract all timeline data
timeline_data <- extract_msk_timeline(msk_data)

# extract for specific cancer types
timeline_data <- extract_msk_timeline(
    msk_data,
    cancer_types = c("Breast Cancer", "Prostate Cancer")
)
```

**Parameters:**
* `mae`: MultiAssayExperiment object from `download_msk_chord()`
* `cancer_types`: Optional character vector of cancer types to include

**Returns:** An `MSKTimeline` S4 object containing:
* `@data`: Named list of timeline data frames
* `@metadata`: Dataset metadata
* `@cancer_types`: Included cancer types
* `@n_patients`: Number of unique patients

**Timeline Types:**
* **Labs**: Tumor markers (CA 15-3, CA 19-9, CEA, PSA)
* **Treatments**: Systemic therapy, surgery, radiation

### 3. Data Filtering and Subsetting

`filter_timeline()`

filter timeline data by multiple criteria.

```R
# filter by patients
filtered <- filter_timeline(
    timeline_obj, 
    patients = c("P-000000XX", "P-000000YY")
)

# filter by date range (in days)
filtered <- filter_timeline(
    timeline_obj,
    date_range = c(-365, 365) # 2 years
)

# filter by event types
filtered <- filter_timeline(
    timeline_obj,
    event_types = c("TREATMENT", "LAB_RESULT")
)

# combine filters
filtered <- filter_timeline(
    timeline_obj,
    patients = c("P-000000XX", "P-000000YY"),
    date_range = c(-365, 365),
    event_types = c("TREATMENT")
)
```

**Parameters:**
* `timeline_obj`: MSKTimeline object
* `patients`: Characyer vector of patient IDs
* `date_range`: Numeric vector of length 2: c(min_days, max_days)
* `event_types`: Character vector of event types to include

**Returns:** Filtered `MSKTimeline` object.

_________________________________________________________________

`extract_by_cancer_type()`

Extract complete patient data for a specific cancer type.

```R
#extract lung cancer cohort (memory efficient, no genomics)
lung_cohort <- extract_by_cancer_type(
    msk_data, #mae
    cancer_type = "Non-Small Cell Lung Cancer"
)

# include genomic data (requires more RAM memory)
lung_cohort <- extract_by_cancer_type(
    msk_data, #mae
    cancer_type = "Non-Small Cell Lung Cancer",
    include_genomics = TRUE
)
```
___________________________________________________________________

`filter_by_genes()`

Filter INCOMMON data for specific genes of interest.

```R
# filter for common cancer genes
filtered <-- filter_by_genes(
    incommon_data,
    genes = c("TP53", "PIK3CA", "KRAS", "PTEN")
)

# check what was retained
print(filtered)
```

**Parameters:**
* `incommon_data`: incommon_ready object
* `genes`: Character vector of gene symbols

**Returns:** Filtered `incommon_ready` object.



___________________________________________________________________

**Parameters:**
* `mae`: MultiAssayExperiment object
* `cancer_type`: Cancer type name (exact match required)
* `include_genomics`: Whether to extract genomic data (default: FALSE)

**Available Cancer Types:**
* `Non-Small Cell Lung Cancer`
* `Breast Cancer`
* `Colorectal Cancer`
* `Prostate Cancer`
* `Pancreatic Cancer`

**Returns:** An `msk_cancer_cohort` object containing:
* `$cancer_type`: Cancer type string
* `$n_patients`: Number of patients
* `$patient_ids`: Patient identifiers
* `$clinical`: Clinical data frame
* `$timeline`: Filtered MSKTimeline object
* `$sample_map`: Sample-to-patient mapping
* `$mae_subset`: Genomic data (if `include_genomics=TRUE`)

__________________________________________________________________

### 4. INCOMMON-Ready Data Preparation
`prepare_incommon_data()`

Prepare data for INCOMMON Bayesian Inference of mutation copy number and multiplicity.

```R
# prepare complete dataset
incommon_data <- prepare_incommon_data(
    cohort = lung_cohort, # msk_cancer_cohort object
    mae = msk_data, 
    max_samples = 100,
    min_purity = 0.1
)

# verify data completeness
verify_incommon_data(incommon_data)
```

**Parameters:**
* `cohort`: msk_cancer_cohort object from `extract_by_cancer_type()`
* `mae`: Original MultiAssayExperiment
* `max_samples`: Maximum samples to extract (default: 100)
* `min_purity`: Minimum tumor purity threshold (default: 0.1)

**Returns:** An `incommon_ready` object containing:
* `$mutations`: Data frame with mutation data
    * `Hugo_symbol`: Gene Symbol
    * `SAMPLE_ID`, `PATIENT_ID`: Identifiers
    * `t_alt_count`: Variant read count
    * `t_ref_count`: Reference read count
    * `t_depth`: Total depth (calculated as `t_alt_count + t_ref_count`)
    * `VAF`: Variant allele frequency (calculated as `t_alt_count / t_depth`)
    * `PURITY`: Tumor purity (fraction)
* `$cna`: Copy number alteration data
* `$purity`: Purity information per sample
* `$parameters`: Dataset parameters

**Data Quality**:
* Automatically calculates depths from alt + ref counts
* filters for complete cases
* Ensures non-negative depths and valid VAF ranges

___________________________________________________________________

`add_copy_number_to_mutations()`

Integrate copy number data into mutation data for INCOMMON analysis.

```R
# add copy number information
incommon_data <- add_copy_number_to_mutations(incommon_data)

# This adds columns:
# - copy_number: log2 ratio from CNA data
# - absolute_cn: Absolute copy number (2 * 2^log2_ratio)
# - k_total: Total copy number for INCOMMON
```

**Parameters:**
* `incommon_data`: Object from `prepare_incommon_data()`

**Returns:** Updated `incommon_ready` object with copy number columns.

___________________________________________________________________

`estimate_eta()`

Estimate the read rate per chromosome copy ($\eta$) for INCOMMON.

```R
# estimate eta from diploid regions
eta <- estimate_eta(incommon_data$mutations)

# use in INCOMMON analysis
# This parameter represents reads per chromosome copy
```

**Parameters:**
* `mutations`: Data frame with mutation data including `k_total`, `VAF`, `t_depth`, `PURITY`

**Returns:** Numeric value of estimated $\eta$.

**Method:**
* Uses diploid regions (k = 2) with moderate VAF (0.2 - 0.6)
* Formula: $\eta = \frac{mean (depth)}{2}$ for diploid regions
* Falls back to overall mean if insufficient diploid mutations

___________________________________________________________________

`check_incommon_readiness()`

Verify all required parameters for INCOMMON analysis are present.

```R
# check readiness
status <- check_incommon_readiness(incommon_data, eta = eta)

# returns information about:
# - Data completeness
# - Copy number availability
# - Parameter ranges
# - Overall readiness status
```

**Parameters:**
* `incommon_data`: Object from `prepare_incommon_data()`
* `eta`: Optional eta value (will note if missing)

**Returns:** Invisible list with readiness status.

____

## Compatibility analysis chordR <-> INCOMMON

### Input format for INCOMMON:

It is needed two data frames:

**genmic_data (mandatory)**:

* sample (Sample ID)
* chr (chromosome)
* from (start position)
* to (end position)
* ref (reference allele)
* alt (alternate allele)
* DP (t_depth)
* NV (variant read count)
* VAF (Variant Allele Frequency)
* gene (name of the gene)
* HGVSp_Short (protein) (optional)

**clinical_data (mandatory)**:

* sample (Sample ID - match with genomic_data)
* purity (tumor purity)

**Optional but helpful**:

* tumor_type (cancer type - ONCOTREE code)
* OS_STATUS (Overall survival status)
* OS_MONTHS (Overall survival time in months)
* SAMPLE_TYPE (Primary/Metastasis)
* MET_COUNT (Number of metastases)
* METASTATIC_SITE (Site of metastasis)

### Data extracted from chordR

From chordR, `prepare_incommon_data()` provides:

`incommon_data$mutations:`

* Hugo_Symbol -> gene
* SAMPLE_ID -> sample
* PATIENT_ID -> (extra, but useful for mapping)
* t_alt_count -> NV
* t_ref_count -> (used to calculate DP/t_depth)
* t_depth -> DP
* VAF -> VAF
* PURITY -> purity (in clinical data)
* seqnames -> chr
* start -> from
* end -> to (if available)
* Reference_Allele -> ref
* Tumor_Seq_allele2 -> alt
* k_total -> (extra, not required)
* HGVSp_Short -> HGVSp_Short

### Install INCOMMON package from GitHub

```R
# install.packages("devtools") # if not already installed
devtools::install_github("caravagnalab/INCOMMON")

library(INCOMMON)
```

## Conversion function

The function `chordr_to_incommon()` will help convert chordR data to INCOMMON format.

```R
# converts data in INCOMMON format
incommon_input <- chordr_to_incommon(
  incommon_data = incommon_data,
  cohort = lung_cohort
)

# verify the results
print(incommon_input)

# verify which analysis can be done
print(incommon_input$conversion_summary$available_analyses)
```


______________________________________________________


___________________________________________________________________

### 5. Visualization 

`plot_timeline()`

Create visual representations of patient timelines.

```R
# plot timelines for specific patients
plot_timeline(
    timeline_data,
    patients = c("P-000000XX", "P-000000YY")
)

# plot first 10 patients (default)
plot_timeline(
    timeline_data
)
```

**Parameters:**
* `timeline_obj`: MSKTimeline object
* `patients`: Character vector of patient IDs (max 20)

**Returns:** A `ggplot2` object

**Features:**
* Color-coded by event type
* Time axis in days from baseline
* Automatic patient selection if not specified 

___________________________________________________________________

### 6. Utility Functions

`get_patient_ids()`

Extract unique patient IDs from a timeline object.

```R
patient_ids <- get_patient_ids(timeline_obj)
```


____________________________________________________________________

## Complete Workflow Examples

### Example 1: Basic Data Exploration

```R
library(chordR)

# download data
msk_data <- download_msk_chord()

# extract timeline data
timeline <- extract_msk_timeline(msk_data)

# get summary
summary(timeline)

# view specific patients
patient_ids <- get_patient_ids(timeline)
head(patient_ids, 10)

# visualize
plot_timeline(timeline, patients = patient_ids[1:5])
```
________________________________________________________________

### Example 2: Cancer-Specific Cohort Analysis

```R
# extract lung cancer cohort
lung_cohort <- extract_by_cancer_type(
    msk_data,
    cancer_type = "Non-Small Cell Lung Cancer"
)

# check cohort size
print(lung_cohort)
# output:
# Cancer Type: Non-Small Cell Lung Cancer
# Patients: 7809
# Clinical variables: 50

# filter timeline to first year post enrollment 
filtered_timeline <- filter_timeline(
    lung_cohort$timeline,
    date_range = c(0, 365)
)

# visualize treatments in first year
plot_timeline(
    filtered_timeline,
    patients = lung_cohort$patient_ids[1:10] # first 10 patients
)
```

________________________________________________________________

### Example 3: INCOMMON Data Preparation

```R
# extract breast cancer with genomics
breast_cohort <- extract_by_cancer_type(
    msk_data,
    cancer_type = "Breast Cancer",
    include_genomics = TRUE
)

# prepare for INCOMMON
incommon_data <- prepare_incommon_data(
    cohort = breast_cohort,
    mae = msk_data,
    max_samples = 50,
    min_purity = 0.2 # higher purity threshold
)

# add copy number data
incommon_data <- add_copy_number_to_mutations(incommon_data)

# estimate eta 
eta <- estimate_eta(incommon_data$mutations)

# verify readiness
check_incommon_readiness(incommon_data, eta = eta)

# filter for genes of interest
cancer_genes <- c("TP53", "PIK3CA", "BRCA1", "BRCA2", "ERBB2")
filtered_data <- filter_by_genes(incommon_data, genes = cancer_genes)

# final verification
verify_incommon_data(filtered_data)

# ready for INCOMMON analysis
```

________________________________________________________________

### Example 4: Multi-Cancer Comparative Analysis

```R
# extract multiple cancer types
cancer_types <- c(
    "Non-Small Cell Lung Cancer",
    "Breast Cancer",
    "Colorectal Cancer"
)

cohorts <- lapply(cancer_types, function(ct) {
    extract_by_cancer_type(msk_data, cancer_type = ct)
})

names(cohorts) <- c("Lung", "Breast", "CRC")

# compare cohort sizes
sapply(cohorts, function(x) x$n_patients)

# combine timelines for comparison
combined_timeline <- do.call(rbind, lapply(cohorts, function(x) {
    x$timeline@data$treatments
}))

# further analysis...
```
________________________________________________________________

## Understanding INCOMMON Integration

### What is INCOMMON?

INCOMMON (**In**tegrated **Co**py number and **M**utation **M**ultiplicity for **On**cology) is a Bayesian method for inferring:

1. **Total copy number (k)** at mutant loci
2. **Mutation multiplicity (m)** - number of copies carrying the mutation
3. **Gene mutant dosage** - fraction of alleles with the mutation

#### Bayesian inference of copy number and mutation multiplicity using INCOMMON

INCOMMON is a Bayesian framework to infer mutation multiplicity and copy number from read counts data. The data consist of a set of $n$ mutations $X = {x_1, ..., x_n}$ from a tumor sample, with each mutation $x_i = (r_i, d_i)$ characterized by the number of reads with the variant $r_i \in \mathbb{Z}^{+}$ and the total number of reads $d_i \in \mathbb{Z}^{+}$ (sequencing depth). INCOMMON jointly models the purity $0 < \pi \leq 1$ and the rate of reads per chromosome copy ($\eta \in \mathbb{R}^{+}$) of the sample, and the total copy number ($k_i \in \mathbb{Z}^{+}$) and mutation multiplicity ($m_i \in \mathbb{Z}^{+}, m_i \leq k_i$ of each mutation in the sample. The joint likelihood of the model is:

$$p(X | \Theta) = \prod_{x_i \in X} p(x_i | \Theta) = \prod_{x_i \in X} p(d_i | \Theta) p(r_i | d_i, \Theta)$$

where $\Theta$ is the vector of model parameters $\Theta = (\pi, \eta, \{k, m\}_{i=1}^{n})$.

In INCOMMON, the total number of reads follows a Poisson distribution:

$$p(d_i | \pi, \eta, k_i) = Poisson(d_i | \lambda)$$

where $\lambda = (1- \pi) 2 \eta + \pi \eta k_i$.

The expected value $\lambda$ is contributed by reads from the normal cells, present in the sample in fraction $1 - \pi$, and from tumour cells present in a fraction $\pi$. The expected number of reads from normal cells, assumed to be diploid, is $2\eta$ across the whole genome. The expected number of reads from tumour cells is $k\eta$, assumed to be proportional to the unknown total copy number.

Given the total read count $d_i$, the number of reads with the variant $r_i$ follows a Binomial distribution:

$$p(r_i | \pi, k_i, m_i) = Poisson(r_i | d_i, \varphi)$$

where $$\varphi = \frac{m \pi}{2(1-\pi) + k \pi}$$

The probability of collecting a read with the variant (success probability) $\varphi$ is equal to the variant allele frequency (VAF) expected for a mutation with copy number configuration ($k, m$), adjusted for the sample purity. INCOMMON embeds prior knowledge of the model parameters $\Theta$ using a prior distribution that can be factorised into:

$$p(\Theta) = p(\pi) p(\eta) p(k, m)$$

INCOMMON uses Markov Chain Monte Carlo (MCMC) sampling to estimate the model posterior distribution $p(\Theta | X)$, that would otherwise be intractable analytically. 


### Why use chordR with INCOMMON?

The MSK-CHORD dataset is ideal for INCOMMON because it provides:
* **Variant read counts** $r_i$: Number of reads supporting the variant
* **Total read depth** $d_i$: Total coverage at the position
* **Tumor purity** $\pi$: From pathology estimates
* **Copy number data**: Log2 ratios from MSK-IMPACT

### INCOMMON Parameters from chordR

After preparing data with `chordR`, we have:

```R
# for each mutation i:
r_i <- incommon_data$mutations$t_alt_count # Variant reads
d_i <- incommon_data$mutations$t_depth     # Total depth
pi  <- incommon_data$mutations$PURITY/100  # Tumor purity (as fraction)
k_i <- incommon_data$mutations$k_total     # Copy number

# sample-level parameter:
n <- estimate_eta(incommon_data$mutations) # Reads per chromosome
```

### Clinical Relevance

Gene mutant dosage has been shown to predict:
* **Overall survival** in multiple cancer types
* **Metastatic propensity**
* **Organ-specific metastatic tropism**
* **Treatment response**

For example, high KRAS mutant dosage in pancreatic cancer is associated with significantly worse prognosis (HR = 3.19, p< 0.0001).

_________________________________________________________________

## Data Structure Reference

### MSKTimeline Object

```R
Class: MSKTimeline
Slots:
    @data           : list of data frames
    @metadata       : list (creation date, source)
    @cancer_types   : character vector 
    @n_patients     : integer

Methods:
    show(object)    # print summary
    summary(object) # detailed summary
```

### msk_cancer_cohort Object

```R
Class: list with class "msk_cancer_cohort"
Elements:
    $cancer_type   : character
    $n_patients    : integer
    $patient_ids   : character vector
    $clinical      : data frame
    $timeline      : MSKTimeline object
    $sample_map    : data frame
    $mae_subset    : MultiAssayExperiment (optional)
```

### incommon_ready Object

```R
Class: list with class "incommon_ready"
Elements:
    $mutations     : data frame with columns:
                     - Hugo_symbol, SAMPLE_ID, PATIENT_ID
                     - t_alt_count, t_ref_count, t_depth
                     - VAF, PURITY
                     - copy_number, absolute_cn, k_total (after add_copy_number)
    $cna           : data frame (gene x sample matrix)
    $purity        : data frame (PATIENT_ID, SAMPLE_ID, PURITY)
    $parameters    : list (n_samples, n_patients, n_mutations, min_purity)
```
_________________________________________________________________

## Tips and Best Practices

### Memory Management
1. Start without genomics: Use `include_genomics = FALSE` for exploration
2. Filter early: Extract specific cancer types before genomic operations
3. Limit samples: Use `max_samples` parameter for INCOMMON preparation
4. Close unused objects: Remove large MAE objects after extraction

```R
# good practice
lung <- extract_by_cancer_type(msk_data, "Non-Small Cell Lung Cancer")
rm(msk_data) # free memory
gc() # garbage collection
```

### Data Quality
1. Verify purity: INCOMMON requires reliable purity estimates
2. Check depth: Ensure adequate sequencing depth (typically > 50x)
3. Filter variants: Remove low-quality or low-VAF mutations
4. Validate results: Always run `verify_incommon_data()` before analysis

**Performance Optimization:**

```R
# for large cohorts, process in batches
batch_size <- 100
n_samples <- nrow(cohort$sample_map)
n_batches <- ceiling(n_samples / batch_size)

results <- list()
for (i in 1:n_batches) {
    # process batch
    batch_data <- prepare_incommon_data(
        cohort, mae, 
        max_samples = batch_size,
    )
    results[[i]] <- analyze_batch(batch_data)
}
```

_________________________________________________________________

## Troubleshooting

### Common Issues

Issue: "Could not download MSK-CHORD data"

```R
# solution: check internet connection and try force refresh
msk_data <- download_msk_chord(force_refresh = TRUE)
```

Issue: "No genomic data available"

```R
# solution: Re-extract with genomics enabled
cohort <- extract_by_cancer_type(
    msk_data, 
    cancer_type = "Breast Cancer", 
    include_genomics = TRUE
)
```

Issue: "Missing required columns"

```R
# solution: Ensure data was prepared correctly
incommon_data <- prepare_incommon_data(cohort, mae)
incommon_data <- add_copy_number_to_mutations(incommon_data)
verify_incommon_data(incommon_data)
```

Issue: "Maximum upload size exceeded" (Shiny App)

```R
# solution: filter data before loading
filtered <- filter_by_genes(incommon_datam, genes = c("TP53", "KRAS"))
# use filtered data in app
```

_________________________________________________________________

## Citation

If you use chordR in your research or MSK-CHORD dataset, please cite: 

**MSK-CHORD Dataset:** 
```
Jee, J., Fong, C., Pichotta, K., et al. (2024). Automated real-world data integration improves cancer outcome prediction. Nature, 636, 728-736. https://doi.org/10.1038/s41586-024-08167-5
```

**INCOMMON Method:**
```
Calonaci, N., Krasniqi, E., Colic, D., et al. (2025). Gene mutant dosage determines prognosis and metastatic tropism in 60,000 clinical cancer samples. Nature Communications (in press).
```

#### Copyright and contacts

Cancer Data Science (CDS) Laboratory.
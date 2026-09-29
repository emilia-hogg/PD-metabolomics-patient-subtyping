# Data Setup

## Required input files

The analysis pipeline requires two external input files that are not stored in this repository:

1. `imputation_2_data`
2. `Full_Patient_sample_information.rds`

Both files must be obtained before running the analysis.

## Metabolomics input

The metabolomics input is: `imputation_2_data`

This is the preprocessed and imputed metabolomics dataset from Natacha's project.

It was generated from the pre-imputation dataset: `log_transformed_preimpute_data_clean_impute_2.rds`
using: `Auto_imputation_2_script.R`

The analysis pipeline uses the resulting `imputation_2_data` file directly, so the imputation does not need to be rerun before running this analysis.

The expected `imputation_2_data` object is a data frame containing:
- 1,781 samples
- 1,122 metabolites

## Clinical input

The clinical and sample information is provided by:

`Full_Patient_sample_information.rds`

The expected clinical dataset contains 1,902 records.

Clinical and metabolomics samples are linked using `Anonymised_sampleID`.

All metabolomics samples should have a corresponding record in the clinical dataset.

## Configure the analysis

Open:

`00_config/config.R`

### 1. Set the project root

Change `project_root` to the location of the repository on the system being used.

For example:

```r
project_root <- "/path/to/Full diss pipeline and outputs"
```

### 2. Set the external input paths

Within the `files` list, set the paths to the metabolomics and clinical input files:

```r
files <- list(
  imputation_2_data = "/path/to/imputation_2_data",
  clinical_rds = "/path/to/Full_Patient_sample_information.rds",
  analysis_cohort = file.path(paths$analysis_dataset, "analysis_cohort.rds"),
  patient_metadata = file.path(paths$analysis_dataset, "patient_metadata.rds"),
  clinical_processed = file.path(paths$analysis_dataset, "clinical_processed.rds"),
  metabolomics_processed = file.path(paths$analysis_dataset, "metabolomics_processed.rds"),
  sample_ids = file.path(paths$analysis_dataset, "sample_ids.rds")
)
```

Replace the example paths with the locations of the input files on the system being used.

The remaining paths used by the analysis are generated from `project_root`.

## Build the analysis cohort

Run:

`01_build_analysis_dataset/01_build_clean_cohort.R`

This script:

- loads `imputation_2_data`;
- loads `Full_Patient_sample_information.rds`;
- matches metabolomics samples to the clinical data using `Anonymised_sampleID`;
- retains patients meeting the study cohort definition;
- requires complete UPDRS-III, MoCA, LEDD, age and sex data;
- calculates disease duration where age at onset is available; and
- saves the analysis cohort and associated quality-control outputs.

A successful run should report:

```text
Final cohort size: 1504
Metabolomics features: 1122
```

## Run the analysis pipeline

After the analysis cohort has been created, run the analysis stages in the following order:

1. `01_build_analysis_dataset`
2. `02_exploratory_analysis`
3. `03_network_construction`
4. `04_snf`
5. `05_sensitivity_analysis`
6. `06_posthoc_analysis`

Within each stage, run scripts in numbered order.

The scripts use the shared paths and analysis settings defined in `00_config/config.R`.

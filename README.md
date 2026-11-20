# CoMMA

**Comparative Modelling of Marker–Trait Associations**

CoMMA is an R workflow for plant breeders working with genotype and phenotype data. It runs genome-wide association studies (GWAS) and genomic prediction, compares prediction accuracy on held-out lines, and uses SHAP to show which genomic regions contribute to predictions. It also handles genotype quality control and imputation, saves fitted models, and creates figures and tables for each trait.

## Quick start in RStudio

The Windows launcher requires **64-bit R 4.6.x** and an internet connection on the first run. It installs missing packages into a library inside the CoMMA folder. [Download R for Windows](https://cran.r-project.org/bin/windows/base/) · [Download RStudio Desktop](https://docs.posit.co/ide/user/download.html)

1. Download this repository with **Code → Download ZIP** and extract it.
2. If you have the TASSEL maize starter ZIP, copy its `input/` folder into the extracted CoMMA folder. It supplies `input/reference.hmp.txt.gz` and `input/phenotypes.csv` for ear height (`EarHT`). To start with your own files, follow [Run your own data](#run-your-own-data) first.
3. Open [Breeding_pipeline.Rproj](Breeding_pipeline.Rproj) in RStudio. Open [Run_pipeline.R](Run_pipeline.R) and click **Source**. The defaults of two preprocessing workers and a 2 GiB memory budget are suitable for the demo.
4. When the console prints `Finished`, open `results/first_run/001_EarHT/report/` for the demo plots and tables.

You can also run `source("Run_pipeline.R")` in the RStudio console with the CoMMA project open.

## Run your own data

Put your files in `input/`. Keep sample IDs consistent across files, including leading zeros, and use the exact trait names from your phenotype table.

| Input | Format |
| --- | --- |
| Genotypes | HapMap (`.hmp.txt` or `.hmp.txt.gz`), or a dosage CSV or TSV with one sample per row. For dosage files, put `sample_id` first, followed by one column per SNP. Use 0, 1 or 2 for the number of copies of the second allele in the marker map and `NA` for missing calls. [Dosage template](templates/reference_calls.tsv) |
| Marker map | Required for dosage files. Supply `marker_id`, `chromosome`, `position_bp` and `alleles` such as `A/G`. The marker IDs must match the genotype columns. HapMap supplies this information in the genotype file. [Marker map template](templates/markers.tsv) |
| Phenotypes | CSV or TSV, either with `sample_id`, `trait` and `value` columns or with `sample_id` and a separate column for each trait. Supply one value per sample and trait. [Phenotype template](templates/phenotypes.tsv) |
| Sample metadata | Optional CSV or TSV with `sample_id` and `source_group` for known cohorts. [Metadata template](templates/sample_metadata.tsv) |

Copy [config/breeding_run.yml](config/breeding_run.yml) to `config/my_run.yml` and adapt it. For example, this runs thousand-kernel weight with two prediction models and all three GWAS methods.

```yaml
reference: input/my_panel.hmp.txt.gz
phenotypes: input/my_phenotypes.csv
output: results/my_run
traits: [TKW]
unit: g

gwas: [MLM, FarmCPU, BLINK]
gwas_candidates: 10
gp: [elastic_net, gaussian_rkhs]
assessment_fraction: 0.25
model_config: config/breeding_pipeline.yml

explain_learners: [gaussian_rkhs]
explain_grouping: coordinate
validate_imputation: true
export_imputed: true
```

For a dosage CSV or TSV, change `reference` and add `markers: input/markers.tsv`. For multiple traits, list them under `traits` or omit that line to run every trait in the table. Use `unit` for a shared plot unit or a mapping such as `unit: {TKW: g, PH: cm}` for different units. CoMMA reads the genotype file once and reuses preprocessing when traits share the same training lines and settings.

In [Run_pipeline.R](Run_pipeline.R), set `run_config <- "config/my_run.yml"` and click **Source**. The input and output paths are relative to the CoMMA folder. Choose a new output folder for each run.

### Choose methods and outputs

| Setting | What it controls |
| --- | --- |
| `gwas` | `MLM`, `FarmCPU` and `BLINK`. Use `gwas: []` to run prediction only. |
| `gp` | `gblup_reml`, `elastic_net`, `random_forest`, `xgboost` and `gaussian_rkhs`. Use `gp: []` for GWAS only. |
| `gwas_candidates` | Maximum size of a follow-up shortlist based on agreement across the three GWAS methods. Omit this setting when using fewer than three GWAS methods. |
| `assessment_fraction` | Fraction of lines held out for prediction assessment. The example holds out 25%. For a GWAS-only run, omit this setting, use `explain_learners: []` and set `validate_imputation: false`. |
| `explain_learners` | Which of the selected GP models to explain with SHAP. Use `[]` to skip SHAP. Choose `explain_grouping: coordinate` for fixed genomic windows or `explain_grouping: ld` to group nearby SNPs based on linkage disequilibrium in the training panel. |
| `validate_imputation` | Mask known genotype calls in held-out lines to estimate imputation accuracy. |
| `export_imputed` | Save the final imputed genotype matrix and its marker map. |

Set `workers` and `memory_gib` at the top of `Run_pipeline.R` to control preprocessing. For a shorter first run on a large panel, select fewer learners or turn off SHAP. SHAP explains the fitted model's predictions. Use the GWAS tables for statistical association results.

## Find your results

Start with `trait_manifest.tsv` in the run folder. It lists the output directory for each trait.

| Within each trait folder | What to look for |
| --- | --- |
| `report/` | Combined figures for the analyses you ran, including Manhattan and QQ plots, PCA, held-out prediction accuracy, imputation accuracy, regional and local SHAP, allele effects, and chromosome plots. |
| `GWAS/` | `GWAS_results.tsv` gives marker results by method. `GWAS_by_marker.tsv` puts the methods side by side. `GWAS_candidate_shortlist.tsv` is available when the shortlist is enabled. |
| `GP/` | The final fitted models, preprocessing and learner settings. Keep the whole folder to predict new lines. |
| `assessment/` | Held-out predictions and accuracy when assessment is enabled. |
| `explanations/` | Regional importance, local contributions and precision checks when SHAP is enabled. |
| `imputed_final/` | `imputed_genotypes.rds` and its matching `markers.tsv` when export is enabled. |

The run folder also contains `hyperparameters.tsv` for GP fits and `first_run_timing.tsv` for the launcher stages. For several traits, `preprocessing_usage.tsv` shows when preprocessing was shared.

## Predict new lines with a saved model

Open the CoMMA project and point to your trait's final `GP/` folder. Your new genotype file must use the same marker IDs and allele coding. If it is a dosage table, pass the matching marker map with `markers = "input/markers.tsv"`. The first line selects the project library created by the Windows launcher if it exists.

```r
if (dir.exists(".breeding-library/R-4.6")) {
  .libPaths(c(".breeding-library/R-4.6", .libPaths()))
}
source("breeding_pipeline.R")

predictions <- predict_breeding_models(
  model = "results/my_run/001_TKW/GP",
  query = "input/new_lines.hmp.txt.gz",
  output = "results/TKW_new_lines"
)
```

## Running from a terminal or an existing R session

On Windows, the PowerShell launcher can use the same run configuration.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\run_pipeline.ps1 -Config config/my_run.yml -Workers 2
```

If you already have a compatible R library, you can source [breeding_pipeline.R](breeding_pipeline.R) from the project folder and call `run_breeding_trait()` or `run_breeding_traits()` directly. For Linux with R 4.5, the [Linux package versions](config/linux_runtime_versions.tsv) list the dependencies used by CoMMA.

## Troubleshooting

- **File missing or sample IDs do not match.** Check the paths in the YAML file and the IDs in the genotype and phenotype files.
- **GAPIT download or checksum error.** The Windows launcher attempts to fetch the pinned GAPIT source. As a fallback, copy `vendor/GAPIT_4.1.0.tar.gz` from the starter ZIP into your CoMMA folder under `vendor/` and retry.
- **Output folder already exists.** Change `output` in your YAML or set `output_folder` in `Run_pipeline.R` to a new path.

See [Third-party notices](THIRD_PARTY_NOTICES.md) for attribution and license details.

# CoMMA

**Comparative Modelling of Marker–Trait Associations**

CoMMA is an R workflow for genome-wide association studies (GWAS), genomic prediction, and SHAP-based interpretation in plant breeding. It supports genotype quality control and imputation, multiple traits, saved prediction models, and combined reports. Open `Breeding_pipeline.Rproj` in RStudio. The repository includes the workflow, settings and input templates.

## First run in RStudio on Windows

Install 64-bit R 4.6.x and RStudio. The optional Windows starter bundle includes the TASSEL EarHT example. To use this source checkout, put your genotype and phenotype files under `input/` and set `config/breeding_run.yml` for your files and traits. Then follow these steps.

1. Open `Breeding_pipeline.Rproj` in RStudio.
2. Open `Run_pipeline.R`. Change `workers` and `memory_gib` at the top if needed. The demo uses two preprocessing workers and a 2 GiB memory budget.
3. Click **Source**. The script sets up the local R library and runs the configured traits. With the starter example, it runs ear height (`EarHT`) on 120 maize lines and 256 SNPs from the public TASSEL tutorial. Follow progress in the RStudio Console.
4. Open the `report/` folder for each trait under the configured output directory. For the starter example, it is `results/first_run/001_EarHT/report/`.

The first run needs internet access. For your own data, copy `config/breeding_run.yml` to a new YAML file and set the genotype and phenotype paths, trait names, and methods there. Omit `sample_metadata` if you have no metadata file. Set `run_config` at the top of `Run_pipeline.R` to that YAML file. Set `output_folder` to a new project-relative folder if you want to rerun with the same config and keep earlier results.

## PowerShell on Windows

Open PowerShell in the extracted starter folder and run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\run_pipeline.ps1 -Workers 2
```

The launcher creates a local package library, installs missing Windows binaries from CRAN and Bioconductor and installs the pinned GAPIT 4.1.0 source after its dependencies. The first launch needs internet access. The Windows starter also contains the small TASSEL EarHT example and the verified GAPIT source archive. The production source archive provides the same launcher. Add your own `input/` files when using that archive.

The command runs genotype preparation, selected GWAS and GP models, held-out prediction, optional regional SHAP and masked-call imputation validation, imputed dosage export, and a combined report. The default demo writes to `results/first_run/001_EarHT/report/`. To rerun or use your own data, edit `config/breeding_run.yml` and choose a new output directory. `-Config config/my_run.yml` selects another run file. Set `markers` in the run file for a dosage CSV or TSV. The exploratory conditional-shuffle analysis uses its own explicit strata and follows the commands below.

`-Workers` controls the maximum number of simultaneous preprocessing workers. The default is one and the workflow reduces the actual count when a phase has fewer tasks or the memory budget allows fewer workers. GP fitting, GWAS and regional explanation still run sequentially in this launcher. The default memory budget is 2 GiB. For example, `-Workers 4 -MemoryGiB 8` requests up to four preprocessing workers and an 8 GiB working-memory budget. Keep that budget within the memory you can spare. If R is installed in a custom location, `-RExecutable` accepts the path to its `R.exe`.

## R library

The core workflow was checked with R 4.5.0 on Linux and R 4.6.1 on Windows. The shortlist and file input releases were checked on CIMMYT TKW in R 4.6.1 on Windows. The Linux package pins are in `config/linux_runtime_versions.tsv`. The earlier public example ran all five GP learners and all three GAPIT methods. A fresh R session reproduced the saved predictions and the Elastic Net, RF and RKHS regional explanations on each platform.

Use the existing validated project library if it is already installed. On Linux, the earlier compatibility setup records its path in `.renv/breeding_library.txt` in the original project. Start R with that library available, then open this extracted folder. For a fresh Linux environment, use the separate R 4.5 compatibility setup bundle, `breeding_pipeline_R45_smoke_20261031_120905.zip`. In an existing Windows project session, load the pinned project library before switching folders. The workflow expects its dependencies to be installed before use.

## Prepare your inputs

Copy the `input/` folder from a CIMMYT TKW input download into the extracted workflow, or make your own `input/` folder using the headers under `templates/`.

| File | What to supply |
| --- | --- |
| HapMap genotype file (`.hmp.txt` or `.hmp.txt.gz`) | Standard tab-separated HapMap with one SNP per row. The workflow reads the map and allele coding from the file. |
| Dosage genotype file (`.csv`, `.tsv` or tab-separated `.txt`) | One line per row. Start with `sample_id` or `Taxa`, then use a column per SNP. Enter 0, 1 or 2. Use `NA` or an empty field for missing calls. |
| Marker map for dosage files | One row per SNP with `marker_id`, `chromosome`, `position_bp` and `alleles`, such as `A/G`. `SNP` or `Name`, `Chromosome`, and `Position` are accepted, as are separate `ref` and `alt` columns. Dosage counts the second allele. |
| Phenotype file | Use `sample_id` or `Taxa` and either `trait` plus `value` columns or one column per trait. The `trait` argument identifies a single trait. The `traits` argument selects traits for a batch. |
| Optional sample metadata file | Use `sample_id` and `source_group` for recorded cohorts. Include the IDs used for fitting and prediction. |

Phenotype, metadata and marker files can be CSV or tab separated. Gzip compression is accepted for the file formats above. Keep IDs as text, including leading zeroes, and use the same spelling in every file. Use physical positions in base pairs, and keep the same marker and allele definitions for training and query inputs. You can also pass the original numeric R matrix and data frames directly.

```r
source("breeding_pipeline.R")

fit <- run_breeding_trait(
  reference = "input/reference.hmp.txt.gz",
  phenotypes = "input/phenotypes.csv", trait = "TKW",
  sample_metadata = "input/sample_metadata.tsv",
  output = "results/TKW_fit", gp = "elastic_net",
  gwas = c("MLM", "FarmCPU", "BLINK"), gwas_candidates = 10L
)

prediction <- predict_breeding_models(
  "results/TKW_fit/GP", query = "input/query.hmp.txt.gz",
  sample_metadata = "input/sample_metadata.tsv",
  output = "results/TKW_prediction"
)

explanation <- explain_breeding_models(
  "results/TKW_fit/GP", query = "input/query.hmp.txt.gz",
  learners = "elastic_net", output = "results/TKW_explanation",
  query_rows = 3L, background_rows = 4L
)
```

For the dosage CSV version of the starter, set `reference = "input/reference_calls.csv"` and `markers = "input/markers.tsv"`. For prediction or explanation, use `query = "input/query_calls.csv"` with the same map. To inspect genotypes before fitting, use `read_breeding_genotypes()` and optionally pass `sample_ids` to load selected lines. Change `TKW` to the exact trait name in your phenotype table. Use a new output folder for each call.

The 96-line starter is intended for a short local run. The full TKW input and its ready-to-run preset are described below. The starter supplies a separate query file for prediction and explanation.

| Analysis | Choices |
| --- | --- |
| GP | `gblup_reml`, `elastic_net`, `random_forest`, `xgboost`, `gaussian_rkhs` |
| GWAS | `MLM`, `FarmCPU`, `BLINK` |
| Regional explanation | Any selected learner saved in the GP folder |

Use an empty vector for either GP or GWAS when you only want the other analysis. Bayesian learners and neural networks are research options outside this source release. The example uses fixed settings. Choose and validate your settings on the appropriate training partitions before fitting a model for deployment.

## Run several traits from one genotype file

Use `run_breeding_traits()` with one phenotype table and one genotype panel. Leave `traits` out to run all traits, or supply their exact names. In a wide phenotype table, keep annotation fields in the sample metadata file.

```r
source("breeding_pipeline.R")

batch <- run_breeding_traits(
  reference = "input/reference.hmp.txt.gz",
  phenotypes = "input/phenotypes_all_traits.csv",
  sample_metadata = "input/sample_metadata.tsv",
  output = "results/multi_trait",
  gp = c("gblup_reml", "gaussian_rkhs"),
  gwas = c("MLM", "FarmCPU", "BLINK"),
  assessment_fraction = 0.25
)
```

The batch imports the genotype file once. Traits with identical training genotypes, marker maps and settings share training-line and SNP QC, imputation and scaling. This applies to assessment and final training panels separately. Traits with different phenotype availability use their own available lines and fit preprocessing on their own training panels. Learner fitting, prediction and association are per trait. Traits run sequentially under the same CPU and memory settings. Use a common `config` for the batch. If traits need different settings, `read_breeding_genotypes()` can load the raw data once for several calls using its `genotype` and `markers` outputs.

For a full panel on the same server used for the TKW run, you can copy just its resource settings into the generic configuration and pass `config = cfg, input_memory_gib = 8` to the batch call.

```r
cfg <- yaml::read_yaml("config/breeding_pipeline.yml")
cfg$resources <- yaml::read_yaml("config/breeding_cimmyt_tkw.yml")$resources
```

`trait_manifest.tsv` lists each output folder. After completion, `preprocessing_usage.tsv` records which fits shared preprocessing. Each trait directory contains the same outputs as a single-trait run, and `batch$traits` contains the individual receipts. Use `unit = c(TKW = "g", PH = "cm")` to label plots for those two selected traits.

The batch-level `hyperparameters.tsv` includes the learner settings for each trait and fit stage.

### Timed CIMMYT run for TKW, grain protein and irrigated yield

Use source release v11 or later and the separate CIMMYT three-trait input archive. Copy its `input/` folder to the workflow folder, or pass its input directory as the second argument. The same 3,484-line raw genotype panel supplies TKW and grain protein. Irrigated yield uses the 2,752 lines with measured values on that panel. Protein and yield each have one additional line outside the shared genotype panel in the older trait-specific evaluation. The run uses the same fixed TKW preset for all three traits, a 25% cohort/genetically balanced assessment holdout, all five GP learners and MLM, FarmCPU and BLINK. The original calibration lines stay in training. The Linux preset allows up to eight one-thread preprocessing workers within the ten-CPU allocation.

Use the pinned library and project-local temporary directories as shown in the full TKW instructions below. The core command records the full shell time in the log and writes detailed timing to the results directory. Pass a project-relative configuration file after the results directory to change the allocation.

```bash
set -o pipefail
{ time R --vanilla --slave -f scripts/run_cimmyt_three_traits.R --args core input results/CIMMYT_three_traits; } 2>&1 | tee cimmyt_three_traits_core.log
```

Read `timing.tsv` for overall and per-trait elapsed seconds, shared genotype import, training and query preprocessing, each GP fit and prediction, GWAS by method and plots. Rows ending in `_total` include their component steps. Use `preprocessing_usage.tsv` to identify reused training preparation. TKW and protein have the same sample availability and can share training preparation at both fit stages. Yield has its own fitted training preparation. Use `trait_manifest.tsv` to find each trait's outputs.

The optional phases run on the completed batch and record wall time per trait and per phase in `optional_timing.tsv`. Each explanation directory also has `explanation_timing.tsv` with per-learner calculation time. Omit the last `explain` argument to explain all five learners. The export phase writes a standalone final imputed dosage matrix for each trait.

```bash
R --vanilla --slave -f scripts/run_cimmyt_three_traits.R --args explain input results/CIMMYT_three_traits random_forest,gaussian_rkhs
R --vanilla --slave -f scripts/run_cimmyt_three_traits.R --args imputation input results/CIMMYT_three_traits
R --vanilla --slave -f scripts/run_cimmyt_three_traits.R --args report input results/CIMMYT_three_traits
R --vanilla --slave -f scripts/run_cimmyt_three_traits.R --args export_imputed input results/CIMMYT_three_traits
```

For the exploratory TKW conditional follow-up after this batch, use the existing TKW launcher with the batch input directory and `results/CIMMYT_three_traits/001_TKW` as the results argument to `conditional_prepare` and `conditional`.

## Choose GWAS regions for follow-up

Set `gwas_candidates = 10L` to add `GWAS/GWAS_candidate_shortlist.tsv`. The option uses all three GWAS methods from the same run. Change the number to adjust the cap, or omit the argument to run GWAS with any selected methods. The shortlist ranks 10-Mb windows by the median of the three method ranks. It keeps one lead SNP per window and checks LD between nearby leads using the original genotype calls. The table includes each lead's three original p-values and counts unresolved nearby LD comparisons. These comparisons had too few shared calls or no genotype variation among those calls. The run receipt returns the shortlist file path.

Use the shortlist to prioritise descriptive follow-up. `GWAS/GWAS_results.tsv` contains one row per marker and method. `GWAS/GWAS_by_marker.tsv` places each method's p-values, BH-adjusted p-values, significance flags and available effect statistics side by side, one row per marker. Each method keeps its original thresholds.

Keep `GP/breeding_model.rds` and `GP/completion.rds` together. Prediction reuses the saved preprocessing and model. In the prediction table, `eligible` and `reason` identify lines excluded by genotype quality checks. For explanations, open `regional_importance.tsv` and `explanation_precision.tsv`. A `review` status identifies an explanation precision screen that needs attention. GWAS p-values are reported separately in `GWAS/GWAS_results.tsv`.

Each GP fit writes `GP/hyperparameters.tsv`. The table identifies configuration settings, pipeline constants, values derived from the training panel and GBLUP's REML estimate. For example, it includes the Elastic Net lambda used by the fit and the actual number of markers sampled at each Random Forest split.

## Make the figures

GWAS writes Manhattan and QQ panels, genomic inflation, the first three genotype PCs, allele associations and chromosome pages under `GWAS/plots/`. The saved-model explanation call writes regional SHAP and local waterfalls under its `plots/` folder. Use `plots = FALSE` to save the numerical evidence for later plotting. Figure export uses the project's `ggplot2` and an R build with Cairo support.

Use `plot_breeding_report()` to put the figures together or draw them again from saved evidence. Supply any combination of completed GWAS, saved explanations, held-out predictions and masked-call imputation evidence. The manifest records which figures were available. No models are fitted when rendering the report.

Reports include `GWAS_by_marker.tsv` when GWAS output is available, and `hyperparameters.tsv` for the final and assessment fits when you pass `trait_output`. The report can recover these settings from earlier saved models. To refresh an existing report, pass a new output directory.

```r
assessment <- read.csv("input/holdout_phenotypes.csv", colClasses = c(sample_id = "character"))

mask <- validate_breeding_imputation(
  "results/TKW_fit/GP", query = "input/query.hmp.txt.gz",
  output = "results/TKW_masked_accuracy", maximum_masks = 50000L
)

figures <- plot_breeding_report(
  output = "results/TKW_figures", trait_output = "results/TKW_fit",
  explanation_output = "results/TKW_explanation",
  predictions = "results/TKW_prediction",
  observed = assessment,
  imputation_accuracy = mask$accuracy_file, unit = "g"
)
```

- For prediction accuracy, supply the measured outcomes for an independent assessment panel. Prediction folders and returned prediction objects carry their query IDs, including lines excluded by QC. For standalone prediction tables the report infers the IDs and checks learner and repeat coverage. Pass `assessment_ids` if you have a separate expected ID list to check. For existing out-of-fold predictions, set `assessment_type = "cross_validation"` and include `fold_id` and, if repeated, `repeat_id`. The two lowest pooled held-out RMSE learners get the featured density panels. Set `featured` to override that choice. Fold variation is drawn when fold predictions are supplied.
- For imputation accuracy, `validate_breeding_imputation()` masks observed calls in at least seven eligible independent query lines and scores the fitted LD-kNN imputer alongside reference-mean and reference-mode baselines. The exact-match percentage uses all masked calls, so missing predictions count against it. The summary also reports exact matches among returned calls and coverage. Dosage RMSE measures dosage error among returned calls on the 0–2 dosage scale. You can regenerate the plot from an existing masked-call `accuracy_summary.tsv` or the verification directory from the native imputation search.
- Local waterfalls select low, median and high examples among the explained lines. They use observed trait values when supplied for all explained lines and otherwise use the learners' mean prediction for sample selection. Signed SHAP contributions are in the waterfalls and in the saved attribution tables.
- Allele boxplots show one point for each line with a recorded genotype and phenotype at the plotted marker. Earlier GWAS reports can recover these points from the saved preprocessor and the original phenotype table passed through `observed`. The reconstructed genotype-class counts and summaries must match the saved ones.
- The TKW launcher uses fixed 10 Mb coordinate groups and a shared training background for the regional SHAP report. `explain_breeding_models()` accepts `region_width_bp` for a new explanation run at another resolution. The optional candidate follow-up below writes its own results table.
- On the chromosome pages, gold lines connect the same significant SNP across GWAS methods, and dashed mint lines connect the same top SHAP window across learners. Direction on the GWAS lanes uses a common per-allele estimate from observed calls with the GWAS fixed covariates. Open symbols have exploratory intervals spanning zero. The overview shows chromosomes with Bonferroni hits when any are available. Otherwise it can show chromosomes with top SHAP windows. Plot 05 shows regional importance across the full explanation map.

### Explain saved models with LD-defined groups

In the v18 source release, choose `grouping = "ld"` in `explain_breeding_models()`. Use a new output directory when explaining a previously saved fit.

```r
explanation <- explain_breeding_models(
  "results/TKW_fit/assessment/GP", query = "input/query.hmp.txt.gz",
  learners = c("elastic_net", "random_forest"),
  output = "results/TKW_explanation_ld", grouping = "ld",
  ld_options = list(scope = "chromosome", decay_r2 = 0.20, group_r2 = 0.20)
)
```

For `Run_pipeline.R` or `run_pipeline.ps1`, set the same choice in `config/breeding_run.yml`.

```yaml
explain_grouping: ld
explain_ld: {scope: chromosome, decay_r2: 0.20, group_r2: 0.20}
```

LD uses squared allele-dosage correlations from the saved discovery genotypes, with at least 20 jointly observed calls per pair. The estimator samples up to 128 pairs per distance bin and chromosome, smooths bin means monotonically and interpolates their crossing of `decay_r2`. Three bins with at least 30 usable pairs each are required for a distance estimate. Distances come from the chromosome where supported, then the pooled panel. The labelled fallback is 250 kb when neither estimate is supported. `scope: panel` uses the pooled curve for all chromosomes. When the curve remains above the threshold, the last supported bin provides the search distance and the status records that boundary. The estimator searches up to 50 Mb by default. Override `max_distance_bp`, `fallback_distance_bp`, `min_joint_calls`, `pairs_per_bin` or `min_pairs_per_bin` under `explain_ld` if needed.

Each group starts at an unassigned SNP and recruits nearby unassigned SNPs on its chromosome within the estimated distance and at or above `group_r2` with the anchor. All fitted SNPs remain in the partition. A stricter comparison can set `group_r2: 0.50` while keeping the decay threshold. The 250 kb fallback is the PLINK clumping-distance reference. PLINK's 0.50 clumping cutoff is a reference for the stricter comparison. [PLINK clumping documentation](https://www.cog-genomics.org/plink/2.0/postproc#clump)

Read `ld_decay.tsv` for sampled-pair counts and curves, `ld_distances.tsv` for the distance and its source on each chromosome, and `ld_group_members.tsv` for the actual SNPs in each group. The chromosome plot shows the span of each group. The regional report matches GWAS hits by actual membership. The fixed-coordinate mode remains available as `grouping = "coordinate"` or `explain_grouping: coordinate`.

## Run the full CIMMYT TKW panel

Extract the v5 or later workflow and copy `input/` from the full TKW v2 input archive into the workflow folder. The v2 metadata adds genotype-derived strata, fingerprint groups and assessment eligibility. The 870 original calibration lines stay in training. The fixed seed selects 871 of the remaining 2,614 lines, balanced by trial cohort and genetic stratum. The resulting holdout model fits preprocessing and GP on its 2,613 training lines before genotype QC. The final GP fit and descriptive GWAS use all 3,484 matched lines, subject to QC. This single holdout measures within-panel prediction on the supplied 12,799-SNP starting panel. The earlier repeated-CV benchmark provides fold variation.

Use the pinned R 4.5 library already set up on the server. From the extracted workflow folder, use the following if it sits directly under the existing Linux project.

```bash
IFS= read -r R_LIBS_USER < ../.renv/breeding_library.txt
export R_LIBS_USER R_LIBS="$R_LIBS_USER"
r_run=$(mktemp -d "$PWD/.tkw_r.XXXXXX")
mkdir "$r_run/tmp" "$r_run/user"
export TEMP="$r_run/tmp" TMP="$r_run/tmp" TMPDIR="$r_run/tmp" R_USER="$r_run/user"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
R --vanilla --slave -e 'library(GAPIT); cat(R.version.string, "\nGAPIT", as.character(packageVersion("GAPIT")), "from", find.package("GAPIT"), "\n")' &&
  R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args core
```

Run the library exports again in each new shell. If GAPIT fails to load, first check that `../.renv/breeding_library.txt` exists and points to the installed R 4.5 library. The GAPIT load check prints the original package error, including any missing startup dependency.

The core phase writes `results/TKW_full/assessment/`, `GP/` and `GWAS/`. After it completes, use the same shell for the optional explanation and imputation phases, followed by the combined report. The report uses whichever optional results are available. Each phase protects existing output. To run the workflow in another output folder, pass it as the third argument, for example `--args core input results/TKW_full_second`.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args explain
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args imputation
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args report
```

To regenerate figures from the same saved run in a fresh directory, use the fourth argument to `report`.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args report input results/TKW_full results/TKW_full/report_v8
```

### Run the optional candidate follow-up

Use the saved `assessment/GP` fit and the original genotypes. With source release v10 or later, generate the TKW inputs using the saved assessment split and the original phenotype and metadata files.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args conditional_prepare input results/TKW_full source_group
```

This uses cohort (`source_group`) as an explicit shuffle-stratum choice for an exploratory run. The command assigns roughly 25% of the registered holdout to explanation and the rest to confirmation fitting. It balances the allocation by cohort and available genetic stratum, keeps fingerprint groups together and assigns roles before loading the phenotypes. Specify another metadata column or a comma-separated combination, such as `source_group,genetic_stratum`, to change the shuffle strata. A final argument changes the explanation fraction. Review `input/conditional_followup/role_counts.tsv` to see how many lines each role has in each stratum.

The preparation command writes three input files under `input/conditional_followup/`. You can also prepare them from your own predefined design. Keep the follow-up roles within the registered assessment IDs and assign roles and strata before using confirmation outcomes.

| File | Contents |
| --- | --- |
| `design.tsv` | `sample_id`, `role` and `shuffle_stratum`. Set the role to `confirmation_fit` or `explanation`. For each confirmation-fit line, give a predeclared stratum based on measured conditions under which phenotypes could be shuffled. Explanation rows can leave the stratum empty. |
| `discovery_phenotypes.tsv` | `sample_id`, `trait`, `value` for the saved assessment model's **training** lines. Include every retained training line. |
| `confirmation_phenotypes.tsv` | `sample_id`, `trait`, `value` for the `confirmation_fit` lines only. Explanation lines supply genotypes only. |

From the updated source folder with the pinned R library loaded, run the following for the ten-SNP follow-up with both statistics.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args conditional input results/TKW_full 10 shap
```

Use `association` instead of `shap` to run conditional association alone. An optional final argument sets a new output directory. Supply the paths to your existing `input/` and saved TKW output if they differ from the example. This phase runs three new discovery-only GWAS scans using the saved preparation and fits small candidate-level XGBoost models when SHAP is selected.

For another trait, use the same file structure and call the public R function. The `query` file can contain both follow-up roles and must carry the same map and allele orientation as the discovery fit. Dosage files also take `markers = "input/markers.tsv"`.

```r
source("breeding_pipeline.R")
followup <- run_breeding_conditional_followup(
  model = "results/TKW_fit/assessment/GP",
  query = "input/reference.hmp.txt.gz",
  design = "input/conditional_followup/design.tsv",
  discovery_phenotypes = "input/conditional_followup/discovery_phenotypes.tsv",
  confirmation_phenotypes = "input/conditional_followup/confirmation_phenotypes.tsv",
  sample_metadata = "input/sample_metadata.tsv", trait = "TKW",
  output = "results/TKW_followup", max_candidates = 10L, shap = TRUE
)
```

The runner freezes the discovery shortlist and local LD controls before loading the confirmation outcomes. It uses original homozygous calls. The shuffle blocks combine the supplied strata with calls at up to two local LD controls and the strongest other shortlisted lead, where available. The XGBoost model uses the candidate and local LD controls as predictors. Both tests use the same phenotype shuffles when both are callable. The default is 999 shuffles and one thread. The candidate table includes raw and Bonferroni-adjusted p-values, callable counts and reasons for untestable candidates. All shortlisted SNPs count toward the adjustment. Association is the primary result and SHAP is an optional sensitivity check.

In `conditional_followup/`, open `conditional_results.tsv`, `discovery_shortlist.tsv` and `shuffle_blocks.tsv`. `conditional_completion.rds` marks a completed run and records the seed and versions. The runner labels its real-data results exploratory. In the TKW confounding simulation, cohort strata alone did not control exposure that varied within cohorts.

Source v12 checks the sum of native TreeSHAP values against predictions using a scale-aware tolerance capped at `0.0001` trait units. `shap_replay_error` records the maximum absolute residual and `shap_replay_max_fraction_of_tolerance` records the largest ratio of residual to its allowed tolerance. If a fit exceeds the bound, the error identifies the marker, shuffle and query line. To rerun after a failed v11 follow-up, use v12 with the existing prepared inputs and saved assessment fit, and pass a new follow-up directory as the final argument.

Source v13 renders `plots/09_conditional_followup.png` after the TKW conditional phase. It shows Bonferroni-adjusted p-values and the share of callable confirmation lines with within-block genotype support. Unavailable tests are shown as dashes. For a completed run, render the same figure directly from its results table without rerunning the analysis. The final argument specifies the saved follow-up folder.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args conditional_plot input results/CIMMYT_three_traits/001_TKW results/CIMMYT_three_traits/001_TKW/conditional_followup_v12
```

To include the plot in a refreshed TKW report, supply a new report directory and the saved follow-up directory. The three-trait report automatically includes follow-up results found in a trait's default `conditional_followup/` directory. In R, `plot_breeding_report()` accepts a `conditional_results` path or data frame.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args report input results/CIMMYT_three_traits/001_TKW results/CIMMYT_three_traits/001_TKW/report_v13 results/CIMMYT_three_traits/001_TKW/conditional_followup_v12
```

Source v14 diagnoses limited within-block variation from a completed follow-up and its original raw genotype file. It checks the saved model, calls and original block assignments, then writes `support/candidate_support.tsv` and `support/genotype_counts_by_block.tsv`. No phenotypes are needed. The table shows 0/2 homozygote counts, missing and heterozygous lead calls, and losses due to incomplete or heterozygous control calls. On the same fully callable lines, it compares strata alone, strata plus the other lead, strata plus local LD controls, and the original full blocks. These comparisons are genotype diagnostics. The reported conditional p-values use the original full blocks.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args conditional_support input results/CIMMYT_three_traits/001_TKW results/CIMMYT_three_traits/001_TKW/conditional_followup_v12
```

### Export the fitted imputed marker matrix

`GP/preprocessor.rds` already contains the completed genotype dosages at `$training$genotype`. Export a standalone RDS and the matching marker map and QC tables from that file.

```bash
R --vanilla --slave -f scripts/run_cimmyt_tkw.R --args export_imputed input results/TKW_full final
```

Read `results/TKW_full/imputed/final/imputed_genotypes.rds` with `readRDS()`. Rows carry retained sample IDs, and columns carry the retained analysis-marker IDs in the same order as `markers.tsv`. Values are unscaled dosages of the second allele and may be fractional. Use `assessment` instead of `final` for the separate assessment-training panel, which uses the preprocessing fitted on its own training lines. For another trait, call `export_breeding_imputed("results/my_trait/GP/preprocessor.rds", "results/my_trait/imputed")`. The exporter reads the saved preparation and does not rerun imputation.

### Full TKW settings

| Setting | Full TKW preset |
| --- | --- |
| Server budget | 10 CPUs total, 2 reserved, up to 8 one-thread preprocessing workers. 252 GB reported RAM, 180 GiB memory budget, 8 GiB input import cap. Fit and explanation learners run sequentially in this entry. |
| Genotype QC | 50% line missingness ceiling, marker missingness below 40% and heterozygosity below 5%, predictor MAF above 0.01 and analysis MAF above 0.05. |
| Native imputation | 60 LD neighbours, 20 donor lines, 2% per-line residual-fill ceiling. The optional accuracy phase masks up to 50,000 verification calls in the registered holdout. |
| GWAS | MLM, FarmCPU and BLINK, two genotype PCs, EYT cohort indicators and a ten-region follow-up shortlist. Bonferroni and BH FDR are both 0.05. |
| GP | All five supported learners with the frozen TKW calibration settings in `config/breeding_cimmyt_tkw.yml`. |
| Explanation | All five saved assessment learners, 24 sampled held-out lines, 8 training background lines, 10 Mb windows and 16 forward/reverse pairs for nonlinear learners. The previous Linux queue had per-model medians of 17, 26 and 38 minutes for RKHS, XGBoost and RF at comparable settings. Allow for roughly 80 minutes of sequential explanation calculation as a planning reference, plus loading and preprocessing. |

For your own panel, pass `assessment_fraction` to `run_breeding_trait()` to create a seeded holdout automatically. The split uses `source_group` and any supplied `genetic_stratum`. Optional `fingerprint_group` keeps groups together and `assessment_eligible = FALSE` keeps a line in training. The registry is saved to `assessment/assessment_split.tsv`. The combined report can find these saved predictions from `trait_output` once you supply `observed`.

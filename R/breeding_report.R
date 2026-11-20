breeding_plot_requirements <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install the project's ggplot2 package to render breeding figures.")
  if (!capabilities("cairo")) stop("Breeding figure export requires an R build with Cairo PNG support.")
}

breeding_plot_learner <- function(x) {
  labels <- c(gblup_reml = "GBLUP", elastic_net = "Elastic Net", random_forest = "Random Forest",
    xgboost = "XGBoost", gaussian_rkhs = "Gaussian RKHS", training_mean = "Training mean", cohort_mean = "Cohort mean")
  ifelse(x %in% names(labels), unname(labels[x]), x)
}

breeding_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(
    text = ggplot2::element_text(colour = "#172A3B"), plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(colour = "#526576", size = 9), legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(colour = "#E1E8ED", linewidth = .25),
    strip.text = ggplot2::element_text(face = "bold"), plot.margin = ggplot2::margin(8, 11, 6, 11))
}

breeding_plot_save <- function(output, name, plots, columns, width, height, title, subtitle = "", caption = "",
                               column_widths = rep(1, columns)) {
  file <- file.path(output, paste0(name, ".png"))
  grDevices::png(file, width = width, height = height, units = "in", res = 240, type = "cairo", bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  n <- ceiling(length(plots) / columns)
  layout <- grid::grid.layout(n + 2L, columns, widths = grid::unit(column_widths, "null"),
    heights = grid::unit.c(grid::unit(.83, "in"), grid::unit(rep(1, n), "null"), grid::unit(.58, "in")))
  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = seq_len(columns)))
  grid::grid.text(title, x = grid::unit(.25, "in"), y = .68, just = "left",
                  gp = grid::gpar(col = "#172A3B", fontsize = 19, fontface = "bold"))
  grid::grid.text(subtitle, x = grid::unit(.26, "in"), y = .27, just = "left", gp = grid::gpar(col = "#526576", fontsize = 10))
  grid::popViewport()
  for (i in seq_along(plots)) {
    vp <- grid::viewport(layout.pos.row = 2L + (i - 1L) %/% columns, layout.pos.col = 1L + (i - 1L) %% columns)
    if (inherits(plots[[i]], "grob")) {
      grid::pushViewport(vp)
      grid::grid.draw(plots[[i]])
      grid::popViewport()
    } else print(plots[[i]], newpage = FALSE, vp = vp)
  }
  grid::pushViewport(grid::viewport(layout.pos.row = n + 2L, layout.pos.col = seq_len(columns)))
  grid::grid.lines(x = c(.02, .98), y = c(.92, .92), gp = grid::gpar(col = "#E1E8ED"))
  grid::grid.text(paste(strwrap(caption, width = floor(width * 16)), collapse = "\n"), x = grid::unit(.26, "in"),
                  y = .41, just = "left", gp = grid::gpar(col = "#526576", fontsize = 8.5, lineheight = 1.04))
  grid::popViewport(2L)
  file
}

breeding_plot_genome <- function(markers) {
  chr <- unique(as.character(markers$chromosome))
  chr <- chr[order(suppressWarnings(as.numeric(chr)), chr, method = "radix")]
  lengths <- as.numeric(tapply(markers$position_bp, markers$chromosome, max)[chr])
  gap <- max(1, sum(lengths) * .002)
  offset <- setNames(c(0, head(cumsum(lengths + gap), -1L)), chr)
  list(offset = offset, center = (offset + lengths / 2) / 1e6,
       divider = (head(offset + lengths, -1L) + gap / 2) / 1e6)
}

breeding_gwas_genotype_points <- function(prep, response, marker_ids) {
  raw <- prep$imputer$reference
  stopifnot(identical(rownames(raw), prep$training_ids), identical(names(response), prep$training_ids),
            all(marker_ids %in% colnames(raw)), all(is.finite(response)))
  do.call(rbind, lapply(marker_ids, function(id) {
    dosage <- raw[, id]
    called <- is.finite(dosage)
    stopifnot(all(dosage[called] %in% 0:2))
    data.frame(marker_id = id, sample_id = names(response)[called], dosage = unname(dosage[called]),
               observed = unname(response[called]))
  }))
}

breeding_gwas_genotype_stats <- function(y) {
  q <- if (length(y)) stats::quantile(y, c(0, .25, .5, .75, 1), names = FALSE) else rep(NA_real_, 5L)
  data.frame(n = length(y), minimum = q[1L], q1 = q[2L], median = q[3L], q3 = q[4L], maximum = q[5L],
             mean = if (length(y)) mean(y) else NA_real_)
}

breeding_gwas_plot_evidence <- function(prep, response, metadata, structure, covariates, scans, shortlist = NULL) {
  leads <- if (!is.null(shortlist)) head(shortlist$lead_marker_id, 10L) else {
    ranked <- scans[is.finite(scans$p_value), ]
    ranked <- ranked[order(ranked$p_value, ranked$marker_id, method = "radix"), ]
    head(unique(ranked$marker_id), 10L)
  }
  hits <- scans$marker_id[which(scans$significant_bonferroni)]
  markers <- prep$markers
  ids <- prep$training_ids
  raw <- prep$imputer$reference
  stopifnot(identical(rownames(raw), ids), identical(names(response), ids),
            identical(as.character(metadata$sample_id), ids), all(c(leads, hits) %in% markers$marker_id))
  scores <- structure$scores[match(ids, structure$scores$sample_id), , drop = FALSE]
  stopifnot(identical(as.character(scores$sample_id), ids), is.null(covariates) || identical(as.character(covariates$Taxa), ids))
  scores$cohort <- metadata$source_group
  covariates <- if (is.null(covariates)) data.frame(row.names = ids) else covariates[, -1L, drop = FALSE]
  estimates <- lapply(unique(c(leads, hits)), function(id) {
    marker <- markers[match(id, markers$marker_id), ]
    dose <- raw[, id]
    called <- is.finite(dose)
    fit_data <- data.frame(observed = unname(response[called]), dosage = dose[called], covariates[called, , drop = FALSE])
    effect <- c(NA_real_, NA_real_, NA_real_)
    if (sum(called) >= 3L && length(unique(fit_data$dosage)) > 1L) {
      fit <- stats::lm(observed ~ ., data = fit_data)
      if (fit$df.residual > 0L && is.finite(stats::coef(fit)["dosage"])) {
        effect <- c(stats::coef(fit)["dosage"], stats::confint(fit, "dosage"))
      }
    }
    data.frame(marker_id = id, alleles = marker$alleles, counted_allele = substr(marker$alleles, 3L, 3L),
               n_called = sum(called), difference = effect[1L], lower = effect[2L], upper = effect[3L])
  })
  estimates <- do.call(rbind, estimates)
  ranked <- estimates[match(leads, estimates$marker_id), ]
  ranked <- ranked[order(-abs(ranked$difference), ranked$marker_id, na.last = NA, method = "radix"), ]
  plot_markers <- head(ranked$marker_id, 3L)
  points <- breeding_gwas_genotype_points(prep, response, plot_markers)
  classes <- lapply(plot_markers, function(id) {
    observed <- points[points$marker_id == id, ]
    alleles <- strsplit(markers$alleles[match(id, markers$marker_id)], "/", fixed = TRUE)[[1L]]
    stopifnot(length(alleles) == 2L)
    genotype <- c(paste0(alleles[1L], "/", alleles[1L], " (reference)"), paste(alleles, collapse = "/"), paste(rep(alleles[2L], 2L), collapse = "/"))
    stats <- lapply(0:2, function(g) breeding_gwas_genotype_stats(observed$observed[observed$dosage == g]))
    data.frame(marker_id = id, dosage = 0:2, genotype = genotype, do.call(rbind, stats))
  })
  list(markers = markers, scores = scores, variance = structure$variance, estimates = estimates,
       lead_ids = leads, genotype_classes = do.call(rbind, classes), genotype_points = points,
       lead_source = if (is.null(shortlist)) "Lowest GWAS p-values across selected methods" else "GWAS candidate shortlist",
       covariates = names(covariates))
}

breeding_gwas_restore_genotype_points <- function(evidence, trait_output, fit, observed) {
  classes <- evidence$genotype_classes
  if (!NROW(classes) || is.null(observed)) return(NULL)
  ids <- as.character(evidence$scores$sample_id)
  if (any(!is.finite(observed$value[match(ids, observed$sample_id)]))) return(NULL)
  directory <- file.path(trait_output, if (length(fit$model_file)) "GP" else "preprocessing")
  path <- gbs_project_path(file.path(directory, "preprocessor.rds"), must_exist = FALSE)
  if (!file.exists(path)) return(NULL)
  checksum <- file.path(directory, "preprocessor_sha256.txt")
  if (length(fit$model_file) && (!file.exists(checksum) ||
      !identical(digest::digest(path, file = TRUE, algo = "sha256"), readLines(checksum, warn = FALSE)))) {
    stop("Saved preprocessing checksum does not match.")
  }
  prep <- readRDS(path)
  if (!inherits(prep, "breeding_pipeline_preprocessor") || !setequal(prep$training_ids, ids) ||
      !identical(prep$markers, evidence$markers)) stop("Saved preprocessing differs from the GWAS panel.")
  response <- gbs_training_response(observed, prep$training_ids)
  points <- breeding_gwas_genotype_points(prep, response, unique(classes$marker_id))
  stats <- do.call(rbind, lapply(seq_len(nrow(classes)), function(i) {
    keep <- points$marker_id == classes$marker_id[i] & points$dosage == classes$dosage[i]
    breeding_gwas_genotype_stats(points$observed[keep])
  }))
  if (!isTRUE(all.equal(stats, classes[names(stats)], tolerance = 1e-10, check.attributes = FALSE))) {
    stop("Original GWAS phenotypes and observed genotype calls are required to restore the individual points.")
  }
  points
}

breeding_plot_gwas <- function(scans, summary, genome, output, trait) {
  methods <- unique(scans$method)
  threshold <- -log10(summary$bonferroni_threshold)
  scans <- scans[is.finite(scans$p_value), ]
  scans$logp <- -log10(pmax(scans$p_value, .Machine$double.xmin))
  scans$position <- (genome$offset[scans$chromosome] + scans$position_bp) / 1e6
  scans$band <- factor(match(scans$chromosome, names(genome$offset)) %% 2L)
  upper <- max(c(threshold, scans$logp))
  panels <- lapply(methods, function(method) {
    x <- scans[scans$method == method, ]
    info <- summary[summary$method == method, ]
    hit <- x[which(x$significant_bonferroni), ]
    labels <- head(hit[order(hit$p_value, hit$marker_id, method = "radix"), ], 3L)
    labels$position_label <- labels$position
    if (nrow(labels) > 1L) {
      indices <- order(labels$position)
      spread <- diff(range(x$position)) * .15
      for (i in seq_along(indices)[-1L]) {
        labels$position_label[indices[i]] <- max(labels$position_label[indices[i]], labels$position_label[indices[i - 1L]] + spread)
      }
      labels$position_label <- labels$position_label - max(0, max(labels$position_label) - max(x$position))
    }
    labels$label_y <- pmin(labels$logp + upper * .08, upper * 1.08)
    manhattan <- ggplot2::ggplot(x, ggplot2::aes(position, logp)) +
      ggplot2::geom_vline(xintercept = genome$divider, colour = "#B7C5D0", linewidth = .28, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(info$bonferroni_threshold), colour = "#C66A17", linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(colour = band), size = .8, alpha = .70) +
      ggplot2::geom_point(data = hit, colour = "#C66A17", size = 1.35) +
      ggplot2::geom_segment(data = labels, ggplot2::aes(xend = position_label, yend = label_y), colour = "#6F8191", linewidth = .24) +
      ggplot2::geom_text(data = labels, ggplot2::aes(position_label, label_y, label = marker_id), size = 2.5, vjust = -.35) +
      ggplot2::scale_colour_manual(values = c("#617A92", "#008987"), guide = "none") +
      ggplot2::scale_x_continuous(breaks = genome$center, labels = names(genome$offset)) +
      ggplot2::scale_y_continuous(limits = c(0, upper * 1.17)) + breeding_plot_theme() +
      ggplot2::labs(title = method, subtitle = sprintf("%d Bonferroni-significant SNPs", nrow(hit)), x = "Chromosome", y = expression(-log[10](italic(P)))) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    qq_data <- data.frame(expected = -log10(stats::ppoints(nrow(x))), observed = sort(x$logp, decreasing = TRUE))
    qq <- ggplot2::ggplot(qq_data, ggplot2::aes(expected, observed)) +
      ggplot2::geom_abline(colour = "#8395A4", linetype = "dashed") + ggplot2::geom_point(colour = "#007F80", size = .75, alpha = .65) +
      breeding_plot_theme() + ggplot2::labs(title = "QQ plot", subtitle = sprintf("Genomic inflation %.3f", info$lambda_gc),
        x = expression(Expected~-log[10](italic(P))), y = expression(Observed~-log[10](italic(P))))
    list(manhattan, qq)
  })
  breeding_plot_save(output, "01_GWAS", unlist(panels, recursive = FALSE), 2L, 12, 1.5 + 2.7 * length(methods),
                     paste("Genome-wide association for", trait), "Manhattan plots, QQ plots and genomic inflation",
                     "Orange points pass Bonferroni. Up to three SNPs with the smallest significant p-values are labelled in each scan.", c(3.2, 1.2))
}

breeding_plot_pca <- function(evidence, output, trait) {
  pcs <- intersect(paste0("PC", 1:3), names(evidence$scores))
  if (length(pcs) < 2L) return(character())
  panels <- lapply(utils::combn(pcs, 2L, simplify = FALSE), function(pc) {
    scores <- evidence$scores
    x <- data.frame(first = scores[[pc[1L]]], second = scores[[pc[2L]]], cohort = scores$cohort)
    label <- function(name) sprintf("%s (%.2f%%)", name, 100 * evidence$variance$proportion[match(name, evidence$variance$component)])
    ggplot2::ggplot(x, ggplot2::aes(first, second, colour = cohort)) +
      ggplot2::geom_hline(yintercept = 0, colour = "#DBE5EB") + ggplot2::geom_vline(xintercept = 0, colour = "#DBE5EB") +
      ggplot2::geom_point(size = 1.35, alpha = .58) + breeding_plot_theme() +
      ggplot2::labs(x = label(pc[1L]), y = label(pc[2L]), colour = "Cohort") +
      ggplot2::theme(legend.position = if (length(unique(x$cohort)) > 1L) "bottom" else "none")
  })
  breeding_plot_save(output, "02_PCA", panels, length(panels), 11.5, 5.2, paste("Genetic structure for", trait),
                     "Leading genotype principal components", "Percentages show variance explained by each component in the GWAS panel.")
}

breeding_plot_effects <- function(evidence, output, trait, unit) {
  x <- evidence$estimates[match(evidence$lead_ids, evidence$estimates$marker_id), ]
  x <- x[is.finite(x$difference), ]
  if (!nrow(x)) return(character())
  x$label <- factor(paste(x$marker_id, x$counted_allele, sep = "  +"),
    levels = paste(x$marker_id, x$counted_allele, sep = "  +")[order(abs(x$difference), x$marker_id)])
  forest <- ggplot2::ggplot(x, ggplot2::aes(difference, label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#526576") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = label), colour = "#007F80", linewidth = .8) +
    ggplot2::geom_point(colour = "#007F80", size = 2.8) + breeding_plot_theme() +
    ggplot2::labs(title = "Adjusted estimates for the leads", x = paste0("Difference per counted allele (", unit, ")"), y = NULL)
  classes <- evidence$genotype_classes
  points <- evidence$genotype_points
  has_points <- NROW(points) > 0L
  limits <- if (NROW(classes)) range(c(classes$minimum, classes$maximum), finite = TRUE)
  panels <- lapply(unique(classes$marker_id), function(id) {
    d <- classes[classes$marker_id == id, ]
    d$class <- factor(d$dosage, levels = 0:2)
    labels <- ifelse(d$n > 0L, sprintf("%s\nn = %d | Mean %.2f\nMedian %.2f", d$genotype, d$n, d$mean, d$median),
                     paste0(d$genotype, "\nn = 0"))
    e <- x[x$marker_id == id, ]
    chart <- ggplot2::ggplot(d[d$n > 0L, ], ggplot2::aes(class))
    if (has_points) {
      observed <- points[points$marker_id == id, ]
      stopifnot(!anyDuplicated(observed$sample_id), identical(as.integer(table(factor(observed$dosage, levels = 0:2))), d$n))
      observed$class <- factor(observed$dosage, levels = 0:2)
      chart <- chart + ggplot2::geom_point(data = observed, ggplot2::aes(y = observed), colour = "#007F80", alpha = .24, size = .75,
        position = ggplot2::position_jitter(width = .21, height = 0, seed = 20261111L))
    }
    chart +
      ggplot2::geom_boxplot(ggplot2::aes(ymin = minimum, lower = q1, middle = median, upper = q3, ymax = maximum),
                           stat = "identity", width = .32, fill = "white", alpha = .58, colour = "#172A3B") +
      ggplot2::geom_point(ggplot2::aes(y = mean), shape = 23, fill = "#A44273", colour = "white", size = 3) +
      ggplot2::scale_x_discrete(labels = labels, drop = FALSE) + ggplot2::coord_cartesian(ylim = limits) + breeding_plot_theme() +
      ggplot2::labs(title = id, subtitle = sprintf("%+.2f %s per %s allele", e$difference, unit, e$counted_allele),
                    x = NULL, y = paste("Observed", trait))
  })
  observation_note <- if (has_points) "Each point is an observed line with a recorded genotype." else
    if (length(panels)) "Individual observations were unavailable for this report." else ""
  breeding_plot_save(output, "07_allele_effects", c(panels, list(forest)), 2L, 13, 8.5, paste("Allele associations for", trait),
    evidence$lead_source, paste(observation_note, "Diamonds show means and box lines show medians. Estimates use the GWAS fixed covariates and exploratory 95% intervals. Reference denotes the first source allele."))
}

breeding_imputation_exact_match <- function(summary, hardcall_error) {
  summary$exact_match_pct_called <- ifelse(summary$n_called > 0, 100 * (1 - hardcall_error), NA_real_)
  summary$exact_match_pct_masked <- ifelse(summary$n_called > 0,
    summary$exact_match_pct_called * summary$n_called / summary$n_masked, 0)
  summary
}

breeding_plot_imputation <- function(input, output) {
  if (is.character(input) && dir.exists(gbs_project_path(input))) {
    directory <- gbs_project_path(input)
    files <- list.files(directory, "^verification_repeat[0-9]+_(overall|by_maf|by_genotype)[.]tsv$")
    parts <- c(overall = "overall", by_maf = "MAF_bin", by_genotype = "true_genotype")
    counts <- table(sub("^verification_repeat[0-9]+_", "", sub("[.]tsv$", "", files)))
    repeat_id <- sub("^verification_repeat([0-9]+)_.*$", "\\1", files)
    part <- sub("^verification_repeat[0-9]+_", "", sub("[.]tsv$", "", files))
    if (!length(files) || anyNA(counts[names(parts)]) || any(table(repeat_id, part) != 1L)) stop("Supply a complete masked verification directory.")
    input <- do.call(rbind, lapply(files, function(file) {
      x <- utils::read.delim(file.path(directory, file))
      type <- sub("^verification_repeat[0-9]+_", "", sub("[.]tsv$", "", file))
      value <- switch(type, overall = rep("all", nrow(x)), by_maf = x$maf_bin, by_genotype = x$true_genotype)
      data.frame(method = "Selected imputer", group = parts[[type]], value = value, n_masked = x$n_masked,
        n_called = x$n_imputed, dosage_RMSE = x$dosage_rmse_called, hardcall_error_called = x$hardcall_error_called)
    }))
  } else if (is.character(input)) input <- utils::read.delim(gbs_project_path(input), colClasses = c(method = "character", value = "character"))
  fields <- c("method", "group", "value", "n_masked", "n_called", "dosage_RMSE", "hardcall_error_called")
  if (!is.data.frame(input) || !all(fields %in% names(input))) stop("Supply a masked-call accuracy summary or the verification output directory.")
  input <- input[input$group %in% c("overall", "MAF_bin", "true_genotype", "truth_genotype", "truth_dosage"), fields]
  input$group[input$group %in% c("truth_genotype", "truth_dosage")] <- "true_genotype"
  called_error <- input$hardcall_error_called[input$n_called > 0L]
  if (!all(c("overall", "MAF_bin", "true_genotype") %in% input$group) || anyNA(input[setdiff(fields, fields[6:7])]) ||
      any(input$n_called < 0 | input$n_called > input$n_masked | input$n_masked <= 0) ||
      any(!is.finite(input$dosage_RMSE[input$n_called > 0L])) || any(!is.finite(called_error) | called_error < 0 | called_error > 1)) {
    stop("Masked-call evidence must include valid overall, MAF and true-genotype summaries.")
  }
  keys <- interaction(input$method, input$group, input$value, drop = TRUE)
  pooled <- do.call(rbind, lapply(split(input, keys), function(x) {
    called <- x$n_called > 0L
    n <- sum(x$n_called)
    data.frame(method = x$method[1L], group = x$group[1L], value = x$value[1L], n_called = n,
               n_masked = sum(x$n_masked), dosage_RMSE = if (n) sqrt(stats::weighted.mean(x$dosage_RMSE[called]^2, x$n_called[called])) else NA_real_,
               error = if (n) stats::weighted.mean(x$hardcall_error_called[called], x$n_called[called]) else NA_real_)
  }))
  pooled$coverage <- pooled$n_called / pooled$n_masked
  pooled <- breeding_imputation_exact_match(pooled, pooled$error)
  overall <- pooled[pooled$group == "overall", ]
  x <- pooled[pooled$group == "true_genotype", ]
  genotype <- c(`0` = "First allele homozygote", `1` = "Heterozygote", `2` = "Second allele homozygote")
  x$value <- factor(x$value, levels = names(genotype), labels = unname(genotype))
  if (anyNA(x$value)) stop("True genotype classes must be coded 0, 1 or 2.")
  class_sizes <- vapply(split(x$n_masked, x$value, drop = TRUE), function(n) {
    if (length(unique(n)) == 1L) n[1L] else NA_real_
  }, numeric(1L))
  class_labels <- setNames(ifelse(is.na(class_sizes), names(class_sizes),
    sprintf("%s (n = %s)", names(class_sizes), format(class_sizes, big.mark = ",", trim = TRUE))), names(class_sizes))
  maf <- pooled[pooled$group == "MAF_bin", ]
  maf_bins <- c("below_0.01", "0.01_to_below_0.05", "0.05_or_more")
  if (all(maf$value %in% maf_bins)) maf$value <- factor(maf$value, levels = maf_bins, labels = c("<1%", "1 to <5%", "5% or more"))
  percentage <- function(y) paste0(y, "%")
  panels <- list(
    ggplot2::ggplot(overall, ggplot2::aes(exact_match_pct_masked, method)) +
      ggplot2::geom_col(width = .52, fill = "#007F80") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", exact_match_pct_masked)), hjust = -.17, size = 3.4) +
      ggplot2::scale_x_continuous(breaks = seq(0, 100, 25), labels = percentage, limits = c(0, 115)) +
      ggplot2::labs(title = "Overall exact genotype matches", x = "Correct among all masked calls", y = NULL) + breeding_plot_theme(),
    ggplot2::ggplot(overall, ggplot2::aes(dosage_RMSE, method)) + ggplot2::geom_point(size = 3, colour = "#007F80", na.rm = TRUE) +
      ggplot2::labs(title = "Overall dosage error", x = "Dosage RMSE", y = NULL) + breeding_plot_theme(),
    ggplot2::ggplot(maf, ggplot2::aes(value, dosage_RMSE, colour = method)) +
      ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(.35), na.rm = TRUE) +
      ggplot2::labs(title = "Dosage error by allele frequency", x = "Minor allele frequency", y = "Dosage RMSE", colour = NULL) + breeding_plot_theme(),
    ggplot2::ggplot(x, ggplot2::aes(exact_match_pct_masked, value, colour = method)) +
      ggplot2::geom_point(size = 2.6, position = ggplot2::position_dodge(.35), na.rm = TRUE) +
      ggplot2::scale_x_continuous(breaks = seq(0, 100, 25), labels = percentage, limits = c(0, 100)) +
      ggplot2::scale_y_discrete(labels = class_labels) +
      ggplot2::labs(title = "Exact genotype matches by true class", x = "Correct among all masked calls in class", y = NULL, colour = NULL) + breeding_plot_theme())
  data.table::fwrite(pooled, file.path(output, "imputation_accuracy_summary.tsv"), sep = "\t")
  breeding_plot_save(output, "03_imputation_accuracy", panels, 2L, 12.5, 7.8,
    "Imputation accuracy on masked observed calls", "Exact genotype matches and allele dosage error",
    sprintf("Exact matches use all masked calls and round continuous dosages to genotype calls. Dosage RMSE uses returned calls. Overall coverage was %.1f-%.1f%%.",
            100 * min(overall$coverage), 100 * max(overall$coverage)))
}

breeding_validate_imputation <- function(model, query, markers = NULL, output, fraction = .05, seed = 1L,
                                         maximum_masks = 50000L, resources = NULL, input_memory_gib = 2, plots = TRUE) {
  if (isTRUE(plots)) breeding_plot_requirements()
  if (!is.numeric(maximum_masks) || length(maximum_masks) != 1L || !is.finite(maximum_masks) ||
      maximum_masks < 1L || maximum_masks > .Machine$integer.max || maximum_masks != floor(maximum_masks)) stop("Supply a positive integer maximum_masks.")
  if (is.character(model)) {
    model <- if (dir.exists(gbs_project_path(model)) || basename(model) == "breeding_model.rds") breeding_load_model(model) else readRDS(gbs_project_path(model))
  }
  prep <- if (inherits(model, "breeding_trait_model")) model$preprocessing else model
  if (!inherits(prep, "breeding_pipeline_preprocessor")) stop("Supply a fitted model or saved preprocessor.")
  if (is.character(query)) {
    imported <- breeding_read_genotypes(query, markers, maximum_memory_gib = input_memory_gib, reference_markers = prep$imputer$markers)
    query <- imported$genotype
    markers <- imported$markers
  } else markers <- breeding_input_markers(markers)
  gbs_validate_calls(query)
  if (length(intersect(rownames(query), prep$training_ids))) stop("Masked validation requires lines outside the imputation reference.")
  gbs_preprocessing_map(markers, prep$predictor_ids, prep$imputer$markers)
  config <- prep$config
  if (!is.null(resources)) config$resources <- resources
  breeding_check_config(config)
  raw <- query[, prep$predictor_ids, drop = FALSE]
  retained <- rowMeans(is.na(raw)) <= config$sample_missing_max
  analysis <- raw[retained, prep$analysis_ids, drop = FALSE]
  if (nrow(analysis) < 7L) stop("Masked verification requires at least seven eligible independent query lines.")
  reference_stats <- prep$imputer$marker_statistics
  plan <- gbs_make_mask_plan(analysis, fraction_per_phase = fraction, seed = seed, maximum_masks_per_phase = maximum_masks,
                             maximum_memory_gib = config$resources$memory_budget_gib, reference_statistics = reference_stats)
  truth <- plan$truth[plan$truth$phase == "verification", ]
  masked <- raw[retained, , drop = FALSE]
  masked[, prep$analysis_ids] <- plan$masked
  if (!dir.create(gbs_project_path(output, must_exist = FALSE), recursive = TRUE, showWarnings = FALSE)) stop("Choose a new imputation validation output folder.")
  grid <- data.frame(ld = config$imputation$ld, k = config$imputation$k, id = "fixed")
  prediction <- predict_gbs_targets_parallel(prep$imputer, masked, truth[c("sample_id", "marker_id")], grid,
                                             config$resources, file.path(output, "checkpoints"))
  idx <- cbind(match(prediction$targets$sample_id, rownames(analysis)), match(prediction$targets$marker_id, colnames(analysis)))
  stat_rows <- match(prediction$targets$marker_id, reference_stats$marker_id)
  stopifnot(!anyNA(idx), !anyNA(stat_rows), !prediction$truth_values_received)
  modes <- vapply(unique(prediction$targets$marker_id), function(id) {
    which.max(tabulate(prep$imputer$reference[, id] + 1L, nbins = 3L)) - 1L
  }, integer(1L))
  values <- list("LD-kNN before residual fill" = prediction$prediction[, 1L],
                 "Reference mean" = 2 * reference_stats$alt_frequency[stat_rows],
                 "Reference mode" = unname(modes[prediction$targets$marker_id]))
  summary <- do.call(rbind, lapply(names(values), function(method) {
    filled <- plan$masked
    filled[idx] <- values[[method]]
    scores <- gbs_score_imputation(filled, truth)
    row <- function(x, group, value) data.frame(method, group, value, n_masked = x$n_masked, n_called = x$n_imputed,
      dosage_RMSE = x$dosage_rmse_called, hardcall_error_called = x$hardcall_error_called, coverage = x$coverage)
    rbind(row(scores$overall, "overall", "all"), row(scores$by_maf, "MAF_bin", scores$by_maf$maf_bin),
          row(scores$by_genotype, "true_genotype", scores$by_genotype$true_genotype))
  }))
  summary <- breeding_imputation_exact_match(summary, summary$hardcall_error_called)
  summary_file <- file.path(output, "accuracy_summary.tsv")
  data.table::fwrite(summary, summary_file, sep = "\t", na = "NA")
  receipt <- list(accuracy_file = summary_file, eligible_lines = rownames(analysis), validation_calls = nrow(truth),
    provenance = plan$provenance, allocation = prediction$allocation)
  gbs_parallel_atomic(receipt, file.path(output, "imputation_validation_completion.rds"))
  if (isTRUE(plots)) receipt$plot_report <- breeding_plot_report(file.path(output, "plots"), imputation_accuracy = summary)
  receipt
}

breeding_prediction_evidence <- function(predictions, ids = NULL) {
  needed <- c("sample_id", "learner", "prediction", "eligible")
  receipt <- NULL
  registered <- NULL
  if (is.character(predictions)) {
    if (dir.exists(gbs_project_path(predictions))) predictions <- file.path(predictions, "predictions.tsv")
    if (basename(predictions) == "predictions.tsv") {
      path <- gbs_project_path(file.path(dirname(predictions), "prediction_receipt.rds"), must_exist = FALSE)
      if (file.exists(path)) {
        receipt <- readRDS(path)
        if (is.null(receipt$assessment_ids)) {
          audit <- breeding_delimited_table(file.path(dirname(predictions), "query_line_QC.tsv"))
          registered <- audit$sample_id
          if (length(registered) != receipt$query_rows) stop("Saved prediction and query QC counts differ.")
        }
      }
    }
    predictions <- breeding_delimited_table(predictions)
    if (!all(needed %in% names(predictions))) stop("Prediction file must contain sample_id, learner, prediction and eligible.")
    predictions$prediction <- breeding_input_numeric(predictions$prediction, "Predictions")
    predictions$eligible <- as.logical(predictions$eligible)
  } else if (is.list(predictions) && !is.data.frame(predictions)) {
    receipt <- predictions$receipt
    predictions <- predictions$predictions
  }
  if (!is.data.frame(predictions) || !all(needed %in% names(predictions)) || !is.numeric(predictions$prediction) ||
      !is.logical(predictions$eligible) || !nrow(predictions)) stop("Supply a prediction table with numeric predictions and logical eligibility.")
  if (!is.null(receipt)) {
    if (!identical(receipt$status, "saved_breeding_models_applied")) stop("Invalid saved prediction receipt.")
    if (!is.null(receipt$assessment_ids)) registered <- receipt$assessment_ids
  }
  if (is.null(ids)) ids <- if (is.null(registered)) sort(unique(as.character(predictions$sample_id))) else registered
  assert_unique_nonmissing(ids, "Assessment sample IDs")
  if (!length(ids) || (!is.null(registered) && !setequal(ids, registered))) stop("Assessment IDs differ from the saved query set.")
  list(predictions = predictions, ids = ids, training_ids = receipt$training_ids,
    source = if (!is.null(registered)) "saved_prediction_query" else "prediction_table_or_explicit_ids")
}

breeding_plot_prediction <- function(predictions, response, ids, type, output, trait, unit, featured, selection_direction) {
  if (type == "holdout" && any(c("fold_id", "repeat_id") %in% names(predictions))) stop("Use assessment_type = 'cross_validation' for fold predictions.")
  if (type == "cross_validation" && !"fold_id" %in% names(predictions)) stop("Cross-validation plots require fold_id.")
  if (!"repeat_id" %in% names(predictions)) predictions$repeat_id <- "1"
  if (!"fold_id" %in% names(predictions)) predictions$fold_id <- "1"
  if (anyNA(predictions[c("repeat_id", "fold_id")])) stop("Fold and repeat IDs cannot be missing.")
  learners <- unique(predictions$learner)
  splits <- split(predictions, predictions$repeat_id)
  score <- lapply(splits, function(x) {
    if (!setequal(x$learner, learners)) stop("Every repeat must include the same learners.")
    mapping <- unique(x[c("sample_id", "fold_id")])
    if (anyDuplicated(mapping$sample_id)) stop("All learners must use the same held-out fold for each line.")
    breeding_score_predictions(x, response, ids, direction = selection_direction)
  })
  d <- predictions[predictions$eligible, ]
  if (!nrow(d)) return(list(files = character(), featured = character()))
  d$observed <- response$value[match(d$sample_id, response$sample_id)]
  d$error_sq <- (d$prediction - d$observed)^2
  by_learner <- split(d, d$learner, drop = TRUE)
  pooled <- do.call(rbind, lapply(names(by_learner), function(id) {
    rows <- by_learner[[id]]
    metrics <- regression_metrics(rows$observed, rows$prediction)
    mean_observed <- mean(rows$observed)
    data.frame(learner = id, rmse = metrics$rmse,
               percent_of_mean = if (mean_observed > 0) 100 * metrics$rmse / mean_observed else NA_real_,
               pearson_r = metrics$pearson_r)
  }))
  top <- do.call(rbind, score)
  pooled$top_decile_recovered <- vapply(pooled$learner, function(id) mean(top$top_fraction_overlap[top$learner == id], na.rm = TRUE), numeric(1L))
  pooled <- pooled[order(pooled$rmse, pooled$learner, method = "radix"), ]
  data.table::fwrite(pooled, file.path(output, "prediction_plot_metrics.tsv"), sep = "\t")
  fitted <- setdiff(pooled$learner, c("training_mean", "cohort_mean"))
  if (!length(fitted)) return(list(files = character(), featured = character()))
  selected <- if (is.null(featured)) head(fitted, 2L) else breeding_release_choices(featured, fitted, "prediction plot learners", allow_empty = FALSE)
  pooled$label <- factor(breeding_plot_learner(pooled$learner), levels = rev(breeding_plot_learner(pooled$learner)))
  folds <- stats::aggregate(error_sq ~ learner + repeat_id + fold_id, d, mean)
  folds$rmse <- sqrt(folds$error_sq)
  folds$label <- factor(breeding_plot_learner(folds$learner), levels = levels(pooled$label))
  density_repeat <- names(splits)[which(names(splits) %in% d$repeat_id)[1L]]
  density <- d[d$repeat_id == density_repeat & d$learner %in% fitted, ]
  n_lines <- length(unique(density$sample_id))
  limits <- range(c(density$observed, density$prediction))
  padding <- if (diff(limits) > 0) diff(limits) * .03 else max(.01, abs(limits[1L]) * .03)
  limits <- limits + c(-padding, padding)
  bins <- max(6L, min(60L, round(sqrt(n_lines / 2.5))))
  breaks <- seq(limits[1L], limits[2L], length.out = bins + 1L)
  density$obs_bin <- findInterval(density$observed, breaks, all.inside = TRUE)
  density$pred_bin <- findInterval(density$prediction, breaks, all.inside = TRUE)
  counts <- stats::aggregate(list(n = rep(1L, nrow(density))), density[c("learner", "obs_bin", "pred_bin")], sum)
  stopifnot(all(tapply(counts$n, counts$learner, sum) == n_lines))
  counts$x <- (breaks[counts$obs_bin] + breaks[counts$obs_bin + 1L]) / 2
  counts$y <- (breaks[counts$pred_bin] + breaks[counts$pred_bin + 1L]) / 2
  count_max <- max(counts$n)
  count_ticks <- unique(pmax(1L, round(pretty(c(1L, count_max), n = 4L))))
  count_ticks <- count_ticks[count_ticks <= count_max]
  plot_density <- function(id) {
    x <- counts[counts$learner == id, ]
    p <- density[density$learner == id, ]
    ggplot2::ggplot(x, ggplot2::aes(x, y, fill = n)) +
      ggplot2::geom_tile(width = diff(breaks)[1L], height = diff(breaks)[1L]) +
      ggplot2::geom_abline(colour = "#C66A17", linetype = "dashed") +
      ggplot2::coord_fixed(xlim = limits, ylim = limits, expand = FALSE) +
      ggplot2::scale_fill_gradient(low = "#DDF0EE", high = "#075C60", limits = c(1, max(2, count_max)), breaks = count_ticks, trans = "sqrt") +
      ggplot2::labs(title = breeding_plot_learner(id), subtitle = sprintf("RMSE %.3g %s", sqrt(mean(p$error_sq)), unit),
                    x = paste("Observed", trait), y = "Held-out prediction", fill = "Lines per bin") + breeding_plot_theme()
  }
  rmse_label <- ifelse(is.finite(pooled$percent_of_mean), sprintf("%.3g / %.1f%%", pooled$rmse, pooled$percent_of_mean), sprintf("%.3g", pooled$rmse))
  correlation_label <- ifelse(is.finite(pooled$pearson_r), sprintf("r = %.2f", pooled$pearson_r), "r unavailable")
  pooled$metric_label <- paste(rmse_label, correlation_label, sep = "  |  ")
  rmse_range <- range(c(pooled$rmse, if (type == "cross_validation") folds$rmse))
  rmse_span <- diff(rmse_range)
  if (rmse_span == 0) rmse_span <- max(rmse_range[2L] * .2, .01)
  metric_column <- rmse_range[2L] + rmse_span * .10
  metric_breaks <- pretty(rmse_range + c(-.05, .10) * rmse_span, n = 4L)
  metric <- ggplot2::ggplot(pooled, ggplot2::aes(rmse, label)) +
    ggplot2::geom_point(data = if (type == "cross_validation") folds else folds[FALSE, ], colour = "#879BAA", size = 1.75,
                        position = ggplot2::position_jitter(height = .12, width = 0, seed = 1)) +
    ggplot2::geom_point(colour = "#007F80", size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = metric_label), x = metric_column, hjust = 0, size = 2.9) +
    ggplot2::scale_x_continuous(limits = rmse_range + c(-.05, .75) * rmse_span,
                              breaks = metric_breaks[metric_breaks <= rmse_range[2L]], expand = c(0, 0)) +
    ggplot2::labs(title = "Pooled RMSE and Pearson r", subtitle = if (type == "cross_validation") "Small points show individual folds" else "Independent assessment set",
                  x = paste0("RMSE (", unit, ")"), y = NULL) + breeding_plot_theme()
  selection <- ggplot2::ggplot(pooled, ggplot2::aes(top_decile_recovered, label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = top_decile_recovered, yend = label), colour = "#D5E3EB", linewidth = .9) +
    ggplot2::geom_point(colour = "#007F80", size = 2.8) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(title = if (selection_direction == "higher") "Top-decile selection overlap" else "Bottom-decile selection overlap",
                  x = "Fraction of the observed selection recovered", y = NULL) + breeding_plot_theme()
  main <- breeding_plot_save(output, "04_prediction_accuracy", c(lapply(selected, plot_density), list(metric, selection)),
    2L, 13, 9.6, paste("Held-out prediction of", trait), sprintf("%s | %d eligible lines in displayed assessment",
      if (type == "holdout") "Independent holdout" else "Cross-validation", n_lines),
    paste(if (type == "cross_validation") "Density panels use the first eligible repeat." else "Density panels use the independent holdout.",
          "RMSE and Pearson r use the same eligible predictions. RMSE percentages divide by the positive observed trait mean."))
  other <- setdiff(fitted, selected)
  if (length(other)) main <- c(main, breeding_plot_save(output, "04b_prediction_other_learners", lapply(other, plot_density),
    min(3L, length(other)), 12.4, 5.2, paste("Held-out predictions for", trait), "Other fitted learners",
    "Axes, bin grid and color scale match the featured prediction plots."))
  list(files = main, featured = selected)
}

breeding_report_regions <- function(directory) {
  file <- gbs_project_path(file.path(directory, "explanation_completion.rds"))
  receipt <- readRDS(file)
  markers <- readRDS(gbs_project_path(file.path(directory, "report_markers.rds")))
  regions <- readRDS(gbs_project_path(file.path(directory, "report_regions.rds")))
  x <- utils::read.delim(gbs_project_path(file.path(directory, "regional_importance.tsv")), colClasses = c(group_id = "character", learner = "character"))
  i <- match(x$group_id, regions$region_id)
  active <- as.character(regions$region_id[regions$active_in_fitted_model])
  if (anyNA(i) || anyDuplicated(x[c("learner", "group_id")]) || !setequal(x$learner, receipt$learners) ||
      !setequal(x$group_id, active) || nrow(x) != length(active) * length(receipt$learners) ||
      any(!is.finite(x$mean_absolute_shapley) | x$mean_absolute_shapley < 0) || any(!is.finite(x$mean_signed_shapley))) {
    stop("Saved regional attribution does not match the explanation receipt and coordinate map.")
  }
  is_ld <- identical(receipt$grouping, "ld")
  membership <- if (is_ld) {
    members <- utils::read.delim(gbs_project_path(file.path(directory, "ld_group_members.tsv")),
      colClasses = c(group_id = "character", marker_id = "character", chromosome = "character", anchor_marker_id = "character"))
    j <- match(members$marker_id, markers$marker_id)
    g <- match(members$group_id, regions$region_id)
    sizes <- table(members$group_id)
    if (anyDuplicated(members$marker_id) || !setequal(members$marker_id, markers$marker_id) ||
        !setequal(members$group_id, active) || anyNA(g) ||
        any(members$chromosome != markers$chromosome[j] | members$position_bp != markers$position_bp[j]) ||
        any(members$chromosome != regions$chromosome[g] | members$anchor_marker_id != regions$anchor_marker_id[g] |
            members$position_bp < regions$region_start[g] | members$position_bp > regions$region_end[g]) ||
        any(x$n_markers != as.integer(sizes[x$group_id]))) stop("LD group membership differs from the saved explanation map.")
    split(members$marker_id, members$group_id)
  }
  fields <- c("chromosome", "region_start", "region_end", if (is_ld) "anchor_marker_id")
  x <- cbind(x, regions[i, fields])
  x <- x[order(match(x$learner, receipt$learners), -x$mean_absolute_shapley, x$group_id, method = "radix"), ]
  x$rank <- as.integer(stats::ave(seq_len(nrow(x)), x$learner, FUN = seq_along))
  list(importance = x, receipt = receipt, markers = markers, directory = directory, membership = membership)
}

breeding_plot_regional <- function(x, hits, has_gwas, genome, output, trait, unit, membership = NULL) {
  x$model <- factor(breeding_plot_learner(x$learner), levels = unique(breeding_plot_learner(x$learner)))
  x$position <- (genome$offset[x$chromosome] + (x$region_start + x$region_end) / 2) / 1e6
  x$band <- factor(match(x$chromosome, names(genome$offset)) %% 2L)
  top <- x[x$rank <= 10L & x$mean_absolute_shapley > 0, ]
  hits <- hits[order(hits$p_value, hits$marker_id, method = "radix"), ]
  is_ld <- !is.null(membership)
  top$gwas_marker <- vapply(seq_len(nrow(top)), function(i) {
    hit <- if (is_ld) which(hits$marker_id %in% membership[[top$group_id[i]]]) else
      which(hits$chromosome == top$chromosome[i] & hits$position_bp >= top$region_start[i] & hits$position_bp <= top$region_end[i])
    if (length(hit)) hits$marker_id[hit[1L]] else NA_character_
  }, character(1L))
  category_name <- if (is_ld) "LD group" else "Window"
  categories <- c(paste(category_name, "contains a Bonferroni SNP"), paste(category_name, "has no Bonferroni SNP"),
    if (is_ld) "Top LD group" else "Top window")
  top$category <- if (has_gwas) ifelse(is.na(top$gwas_marker), categories[2L], categories[1L]) else rep(categories[3L], nrow(top))
  group_labels <- if (is_ld) top$anchor_marker_id else sprintf("%s %.2g-%.2g Mb", top$chromosome, top$region_start / 1e6, top$region_end / 1e6)
  top$label <- ifelse(is.na(top$gwas_marker), group_labels, top$gwas_marker)
  maximum <- max(x$mean_absolute_shapley)
  height <- if (maximum > 0) maximum else 1
  labels <- top[FALSE, ]
  if (nrow(top)) {
    labels <- do.call(rbind, lapply(split(top, list(top$learner, top$category), drop = TRUE), head, 2L))
  }
  labels$x_label <- labels$position
  labels$y_label <- labels$mean_absolute_shapley + height * .12
  if (nrow(labels) > 1L) {
    spread <- diff(range(x$position)) * .12
    for (indices in split(seq_len(nrow(labels)), labels$learner)) {
      indices <- indices[order(labels$position[indices])]
      for (i in seq_along(indices)[-1L]) labels$x_label[indices[i]] <- max(labels$x_label[indices[i]], labels$x_label[indices[i - 1L]] + spread)
      labels$x_label[indices] <- labels$x_label[indices] - max(0, max(labels$x_label[indices]) - max(x$position))
    }
  }
  chart <- ggplot2::ggplot(x, ggplot2::aes(position, mean_absolute_shapley)) +
    ggplot2::geom_vline(xintercept = genome$divider, colour = "#CEDBE4", linewidth = .26, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = band), size = .79, alpha = .64) +
    ggplot2::scale_colour_manual(values = c("#617A92", "#008987"), guide = "none") +
    ggplot2::geom_point(data = top, ggplot2::aes(shape = category, fill = category), size = 2.15, colour = "#172A3B") +
    ggplot2::geom_segment(data = labels, ggplot2::aes(xend = x_label, yend = y_label), colour = "#8395A4", linewidth = .25) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x_label, y_label, label = label), colour = "#172A3B", vjust = -.25, size = 2.65) +
    ggplot2::scale_shape_manual(values = setNames(c(24, 21, 22), categories)) +
    ggplot2::scale_fill_manual(values = setNames(c("#EB9D42", "white", "#6CD0BA"), categories)) +
    ggplot2::facet_grid(model ~ .) + ggplot2::scale_x_continuous(breaks = genome$center, labels = names(genome$offset)) +
    ggplot2::scale_y_continuous(limits = c(0, maximum + height * .28)) + breeding_plot_theme() +
    ggplot2::labs(x = "Chromosome", y = paste0("Mean absolute regional SHAP (", unit, ")"), fill = NULL, shape = NULL)
  file <- breeding_plot_save(output, "05_regional_importance", list(chart), 1L, 12.4, max(4.5, 2.35 + 1.25 * length(unique(x$learner))),
    paste("Regional SHAP importance for", trait, "prediction"), if (is_ld) "LD-defined SHAP contributions" else "Coordinate-window SHAP contributions", if (is_ld)
      "Labels show a significant member SNP or the group anchor. Group position is the midpoint of its marker span."
    else if (has_gwas)
      "Up to two windows in each category are labelled per learner. A shared position is descriptive. SHAP importance ranks predictive contribution."
    else "Up to two top windows are labelled per learner. SHAP importance ranks predictive contribution.")
  data.table::fwrite(top[setdiff(names(top), c("model", "position", "band", "label"))], file.path(output, "top_SHAP_windows.tsv"), sep = "\t")
  list(file = file, top = top)
}

breeding_plot_local <- function(explanation, response, featured, output, trait, unit) {
  learners <- explanation$receipt$learners
  chosen <- if (length(intersect(featured, learners))) head(intersect(featured, learners), 2L) else head(learners, 2L)
  saved <- setNames(lapply(chosen, function(id) {
    result <- readRDS(gbs_project_path(file.path(explanation$directory, id, "explanations.rds")))
    if (!is.matrix(result$contributions) || !identical(rownames(result$contributions), result$query_ids) ||
        !identical(result$query_ids, explanation$receipt$explained_query_ids) ||
        !identical(result$background_ids, explanation$receipt$background_ids) ||
        any(!is.finite(c(result$contributions, result$baseline, result$prediction))) ||
        max(abs(rowSums(result$contributions) + result$baseline - result$prediction)) > 1e-8) stop("Saved local explanations fail identity or additivity checks.")
    result
  }), chosen)
  ids <- explanation$receipt$explained_query_ids
  measured <- if (is.null(response)) rep(NA_real_, length(ids)) else response$value[match(ids, response$sample_id)]
  prediction <- rowMeans(do.call(cbind, lapply(saved, `[[`, "prediction")))
  is_ld <- !is.null(explanation$membership)
  windows <- unique(explanation$importance[c("group_id", "chromosome", "region_start", "region_end", if (is_ld) "anchor_marker_id")])
  stopifnot(!anyDuplicated(windows$group_id))
  window_labels <- if (is_ld) setNames(paste(windows$chromosome, windows$anchor_marker_id), windows$group_id) else
    setNames(sprintf("%s %s-%s Mb", windows$chromosome,
    format((windows$region_start - 1) / 1e6, trim = TRUE, scientific = FALSE),
    format(windows$region_end / 1e6, trim = TRUE, scientific = FALSE)), windows$group_id)
  observed_order <- all(is.finite(measured))
  order_ids <- order(if (observed_order) measured else prediction, ids, method = "radix")
  positions <- if (length(ids) == 1L) c(median = 1L) else if (length(ids) == 2L) c(low = 1L, high = 2L) else
    setNames(round(1 + c(.10, .50, .90) * (length(ids) - 1L)), c("low", "median", "high"))
  plots <- lapply(names(positions), function(example) {
    i <- order_ids[positions[[example]]]
    evidence <- lapply(chosen, function(learner) {
      result <- saved[[learner]]
      value <- result$contributions[i, ]
      keep <- head(order(-abs(value), names(value), method = "radix"), 6L)
      shown <- value[keep]
      labels <- unname(window_labels[names(shown)])
      stopifnot(!anyNA(labels))
      if (length(value) > length(keep)) {
        shown <- c(shown, sum(value[-keep]))
        labels <- c(labels, if (is_ld) "Other LD groups (sum)" else "Other windows (sum)")
      }
      cumulative <- result$baseline[i] + c(0, cumsum(shown))
      stopifnot(abs(tail(cumulative, 1L) - result$prediction[i]) <= 1e-8)
      list(baseline = result$baseline[i], prediction = result$prediction[i], rows = data.frame(
        window = labels, row = rev(seq_along(shown)), contribution = unname(shown), start = head(cumulative, -1L),
        end = tail(cumulative, -1L), direction = ifelse(shown >= 0, "Raises prediction", "Lowers prediction")))
    })
    extent <- range(c(measured[i], unlist(lapply(evidence, function(e) c(e$rows$start, e$rows$end)))), finite = TRUE)
    padding <- if (diff(extent) > 0) .20 * diff(extent) else max(.001, abs(extent[1L]) * .025)
    panels <- Map(function(learner, e) {
      steps <- e$rows
      steps$label <- ifelse(abs(steps$contribution) < .005 & steps$contribution != 0,
        sprintf("%+.2g", steps$contribution), sprintf("%+.2f", steps$contribution))
      chart <- ggplot2::ggplot(steps) +
        ggplot2::geom_vline(xintercept = e$baseline, colour = "#526576", linetype = "dashed", linewidth = .45)
      if (is.finite(measured[i])) chart <- chart + ggplot2::geom_vline(xintercept = measured[i], colour = "#A44273", linetype = "dotted", linewidth = .7)
      chart +
        ggplot2::geom_segment(data = head(steps, -1L), ggplot2::aes(x = end, xend = end, y = row - .29, yend = row - .71),
                             colour = "#526576", linewidth = .38) +
        ggplot2::geom_rect(ggplot2::aes(xmin = pmin(start, end), xmax = pmax(start, end), ymin = row - .29, ymax = row + .29,
                                      fill = direction)) +
        ggplot2::geom_point(data = tail(steps, 1L), ggplot2::aes(end, row), size = 2.25, colour = "#172A3B") +
        ggplot2::geom_text(ggplot2::aes(x = pmax(start, end) + padding * .1, y = row, label = label), hjust = 0, size = 3.05) +
        ggplot2::scale_fill_manual(values = c("Raises prediction" = "#007F80", "Lowers prediction" = "#C66A17"), name = NULL) +
        ggplot2::scale_y_continuous(breaks = steps$row, labels = steps$window) +
        ggplot2::coord_cartesian(xlim = extent + padding * c(-.3, 1)) + breeding_plot_theme() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) +
        ggplot2::labs(title = breeding_plot_learner(learner), subtitle = sprintf("Baseline %.2f | Prediction %.2f%s", e$baseline, e$prediction,
          if (is.finite(measured[i])) sprintf(" | Error %+.2f", e$prediction - measured[i]) else ""),
          x = paste0("Cumulative prediction (", unit, ")"), y = NULL)
    }, chosen, evidence)
    subtitle <- sprintf("%s example | %s | Selected by %s%s", tools::toTitleCase(example), ids[i],
      if (observed_order) "observed trait" else "mean prediction", if (is.finite(measured[i])) sprintf(" | Observed %.3g %s", measured[i], unit) else "")
    file <- breeding_plot_save(output, paste0("06_local_explanation_", example), panels, length(panels), 13, 5.8,
      paste("Local SHAP attribution for", trait, "prediction"), subtitle,
      "Bars show signed contributions in trait units for up to six windows and the sum of the others. The dashed line marks the baseline, the dot marks the fitted prediction and the magenta line marks the observed value when supplied.")
    list(file = file, data = data.frame(example = example, sample_id = ids[i], learner = chosen, observed = measured[i],
                                      prediction = vapply(evidence, `[[`, numeric(1L), "prediction")))
  })
  data.table::fwrite(do.call(rbind, lapply(plots, `[[`, "data")), file.path(output, "local_explanation_selection.tsv"), sep = "\t")
  vapply(plots, `[[`, character(1L), "file")
}

breeding_conditional_plot_data <- function(results) {
  if (is.character(results)) results <- breeding_delimited_table(results)
  numeric_fields <- c("priority", "n_confirmation", "n_informative", "association_p_bonferroni", "shap_p_bonferroni")
  fields <- c(numeric_fields, "lead_marker_id", "association_status", "shap_status")
  if (!is.data.frame(results) || !nrow(results) || !all(fields %in% names(results))) {
    stop("Supply a conditional_results.tsv table with candidate status, adjusted p-values and support counts.")
  }
  for (field in numeric_fields) results[[field]] <- breeding_input_numeric(results[[field]], field)
  assert_unique_nonmissing(results$lead_marker_id, "Conditional candidate SNP IDs")
  if (anyNA(results[numeric_fields[1:3]]) || !setequal(results$priority, seq_len(nrow(results))) ||
      anyDuplicated(results$priority) || any(results$n_confirmation != floor(results$n_confirmation) |
      results$n_informative != floor(results$n_informative)) || any(results$n_confirmation < 0 | results$n_informative < 0 |
      results$n_informative > results$n_confirmation)) stop("Conditional ranks or support counts are inconsistent.")
  for (method in c("association", "shap")) {
    status <- results[[paste0(method, "_status")]]
    p <- results[[paste0(method, "_p_bonferroni")]]
    if (anyNA(status) || any((status == "tested") != !is.na(p)) || any(p <= 0 | p > 1, na.rm = TRUE)) {
      stop("Conditional test status and adjusted p-values disagree for ", method, ".")
    }
  }
  results[order(results$priority), , drop = FALSE]
}

breeding_plot_conditional <- function(results, output, trait) {
  n <- nrow(results)
  tested <- sum(results$association_status == "tested")
  paired <- any(results$shap_status != "disabled")
  methods <- c("Conditional association", if (paired) "Conditional-shuffle SHAP")
  colors <- setNames(c("#A85B0C", "#007E72")[seq_along(methods)], methods)
  results$row <- n + 1L - results$priority
  rows <- do.call(rbind, lapply(seq_along(methods), function(i) {
    status <- if (i == 1L) results$association_status else results$shap_status
    p <- if (i == 1L) results$association_p_bonferroni else results$shap_p_bonferroni
    data.frame(row = results$row + if (paired) c(.15, -.15)[i] else 0, method = methods[i],
      score = -log10(p), p = p, status = status)
  }))
  called <- rows[rows$status == "tested", , drop = FALSE]
  unavailable <- rows[rows$status != "tested", , drop = FALSE]
  max_score <- max(c(2, called$score))
  threshold <- -log10(.05)
  theme <- breeding_plot_theme() + ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(colour = "#E1E8ED"),
    axis.text.y = ggplot2::element_text(colour = "#172A3B"))
  y_axis <- function(labels) ggplot2::scale_y_continuous(breaks = results$row, labels = labels,
    limits = c(.5, n + .5), expand = c(0, 0))
  axis_p <- c(1, .1, .05, .01, .001, .0001, .00001)
  axis_p <- axis_p[-log10(axis_p) <= max_score + .01]
  tests <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = threshold, colour = "#6E7E89", linetype = "dashed", linewidth = .45) +
    ggplot2::geom_point(data = called, ggplot2::aes(score, row, colour = method, shape = method), size = 2.8) +
    ggplot2::geom_text(data = unavailable, ggplot2::aes(y = row), label = "-", x = -.35, colour = "#526576", size = 4) +
    ggplot2::scale_x_continuous(breaks = -log10(axis_p), labels = format(axis_p, scientific = FALSE, trim = TRUE),
      limits = c(-.55, max_score + .25), expand = c(0, 0)) +
    y_axis(paste0(results$priority, "   ", results$lead_marker_id)) +
    ggplot2::scale_colour_manual(values = colors, limits = methods, drop = FALSE) +
    ggplot2::scale_shape_manual(values = setNames(c(16, 17)[seq_along(methods)], methods), limits = methods, drop = FALSE) +
    ggplot2::labs(title = "Conditional evidence", x = "Bonferroni-adjusted p (log scale)", y = NULL, colour = NULL, shape = NULL) + theme

  results$informative_share <- ifelse(results$n_confirmation > 0, results$n_informative / results$n_confirmation, 0)
  results$counts <- paste0(results$n_informative, " / ", results$n_confirmation)
  support <- ggplot2::ggplot(results, ggplot2::aes(y = row)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, yend = row), linewidth = 2.5, colour = "#E1E8ED") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = informative_share, yend = row), linewidth = 2.5, colour = "#337C9B") +
    ggplot2::geom_vline(xintercept = .4, colour = "#6E7E89", linetype = "dashed", linewidth = .45) +
    ggplot2::geom_text(ggplot2::aes(label = counts), x = 1.05, hjust = 0, colour = "#172A3B", size = 3.1) +
    ggplot2::scale_x_continuous(breaks = c(0, .4, 1), labels = c("0%", "40%", "100%"), limits = c(0, 1.36), expand = c(0, 0)) +
    y_axis(rep("", n)) + ggplot2::labs(title = "Support within shuffle blocks", x = "Informative / callable confirmation lines", y = NULL) + theme
  # Use the same legend space so the candidate rows line up across panels.
  legend <- ggplot2::ggplotGrob(tests)
  support <- ggplot2::ggplotGrob(support)
  bottom <- which(legend$layout$name == "guide-box-bottom")
  support_bottom <- which(support$layout$name == "guide-box-bottom")
  if (length(bottom) && length(support_bottom)) {
    support$heights[support$layout$t[support_bottom]] <- legend$heights[legend$layout$t[bottom]]
    support$heights[support$layout$t[support_bottom] - 1L] <- legend$heights[legend$layout$t[bottom] - 1L]
  }
  breeding_plot_save(output, "09_conditional_followup", list(legend, support), 2L, 12.4, max(5.5, 2.5 + .35 * n),
    paste("Conditional SNP follow-up for", trait), sprintf("Exploratory analysis | %d of %d candidates testable for association", tested, n),
    "Dash marks an unavailable test. Dashed guides show adjusted p = 0.05 and the minimum 40% informative share. Adjustment includes every shortlisted SNP.",
    column_widths = c(1.16, 1))
}

breeding_plot_report <- function(output, trait_output = NULL, gwas_output = NULL, explanation_output = NULL,
                                predictions = NULL, observed = NULL, assessment_ids = NULL,
                                assessment_type = c("holdout", "cross_validation"), imputation_accuracy = NULL,
                                trait = NULL, unit = "trait units", featured = NULL, selection_direction = "higher",
                                conditional_results = NULL) {
  breeding_plot_requirements()
  assessment_type <- match.arg(assessment_type)
  if (!is.null(featured) && (!is.character(featured) || length(featured) < 1L || length(featured) > 2L ||
      anyNA(featured) || anyDuplicated(featured))) stop("Choose one or two featured learner IDs.")
  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) stop("Supply one trait unit label.")
  fit_trait <- NULL
  fit <- NULL
  if (!is.null(trait_output)) {
    fit <- readRDS(gbs_project_path(file.path(trait_output, "trait_completion.rds")))
    if (!identical(fit$status, "complete")) stop("Supply a completed trait output folder.")
    fit_trait <- fit$trait
    if (is.null(gwas_output) && length(fit$gwas_file)) gwas_output <- file.path(trait_output, "GWAS")
    if (is.null(predictions) && !is.null(observed) && !is.null(fit$assessment)) {
      predictions <- file.path(trait_output, "assessment", "predictions")
    }
  }
  prediction_input <- if (!is.null(predictions)) breeding_prediction_evidence(predictions, assessment_ids)
  if (!is.null(prediction_input) && assessment_type == "holdout" && length(fit$model_file) && is.null(prediction_input$training_ids)) {
    prediction_input$training_ids <- breeding_delimited_table(file.path(trait_output, "GP", "training_line_QC.tsv"))$sample_id
  }
  if (!is.null(prediction_input$training_ids)) {
    assert_unique_nonmissing(prediction_input$training_ids, "Prediction model training panel IDs")
    if (length(intersect(prediction_input$ids, prediction_input$training_ids))) stop("Holdout IDs overlap the prediction model training panel.")
  }
  association <- evidence <- explanation <- NULL
  if (!is.null(gwas_output)) {
    association <- readRDS(gbs_project_path(file.path(gwas_output, "GWAS_results.rds")))
    path <- gbs_project_path(file.path(gwas_output, "report_data.rds"), must_exist = FALSE)
    if (file.exists(path)) evidence <- readRDS(path)
  }
  if (!is.null(explanation_output)) explanation <- breeding_report_regions(explanation_output)
  labels <- unique(c(trait, fit_trait, association$summary$trait, explanation$receipt$trait))
  if (length(labels) > 1L || anyNA(labels) || any(!nzchar(labels))) stop("Report inputs must belong to one named trait.")
  if (length(labels) == 1L) trait <- labels[[1L]] else trait <- "trait"
  if (!is.null(predictions) && is.null(observed)) stop("Observed assessment phenotypes are required for prediction accuracy plots.")
  if (!is.null(observed)) {
    observed <- breeding_input_phenotypes(observed, trait)
    if (!is.data.frame(observed) || !all(c("sample_id", "value") %in% names(observed)) || !is.numeric(observed$value)) {
      stop("Supply observed phenotypes with sample_id and numeric value.")
    }
    if ("trait" %in% names(observed)) observed <- observed[which(observed$trait == trait), , drop = FALSE]
    assert_unique_nonmissing(as.character(observed$sample_id), "Observed report sample IDs")
  }
  if (!is.null(evidence) && is.null(evidence$genotype_points) && length(fit$gwas_file) &&
      identical(gbs_project_path(gwas_output), gbs_project_path(file.path(trait_output, "GWAS")))) {
    evidence$genotype_points <- breeding_gwas_restore_genotype_points(evidence, trait_output, fit, observed)
  }
  fields <- c("marker_id", "chromosome", "position_bp")
  scans <- association$results
  gwas_markers <- unique(scans[fields])
  maps <- rbind(gwas_markers, explanation$markers[fields])
  if (!is.null(maps) && nrow(maps)) {
    maps$chromosome <- as.character(maps$chromosome)
    maps <- unique(maps)
    if (anyDuplicated(maps$marker_id)) stop("GWAS and explanation marker coordinates differ.")
  }
  conditional <- if (!is.null(conditional_results)) breeding_conditional_plot_data(conditional_results)
  if (!length(fit$model_file) && all(vapply(list(association, explanation, predictions, imputation_accuracy, conditional), is.null, logical(1L)))) {
    stop("Supply report evidence to plot.")
  }
  if (!dir.create(gbs_project_path(output, must_exist = FALSE), recursive = TRUE, showWarnings = FALSE)) stop("Choose a new figure output folder.")
  files <- list(GWAS = character(), PCA = character(), imputation = character(), prediction = character(),
    regional_importance = character(), local_explanation = character(), allele_effects = character(), chromosome_overview = character())
  if (!is.null(conditional)) files$conditional_followup <- breeding_plot_conditional(conditional, output, trait)
  tables <- list()
  if (length(fit$model_file)) {
    tables$hyperparameters <- file.path(output, "hyperparameters.tsv")
    data.table::fwrite(breeding_trait_parameters(trait_output, fit), tables$hyperparameters, sep = "\t")
  }
  hits <- data.frame(marker_id = character(), method = character(), chromosome = character(), position_bp = numeric(), p_value = numeric())
  methods <- if (is.null(scans)) character() else unique(scans$method)
  if (!is.null(scans)) {
    tables$gwas_by_marker <- file.path(output, "GWAS_by_marker.tsv")
    data.table::fwrite(breeding_gwas_by_marker(scans), tables$gwas_by_marker, sep = "\t", na = "NA")
    genome <- breeding_plot_genome(gwas_markers)
    files$GWAS <- breeding_plot_gwas(scans, association$summary, genome, output, trait)
    hits <- scans[which(scans$significant_bonferroni), names(hits)]
  }
  if (!is.null(evidence)) {
    files$PCA <- breeding_plot_pca(evidence, output, trait)
    files$allele_effects <- breeding_plot_effects(evidence, output, trait, unit)
    data.table::fwrite(evidence$estimates, file.path(output, "allele_association_estimates.tsv"), sep = "\t")
    if (!is.null(evidence$genotype_classes)) data.table::fwrite(evidence$genotype_classes, file.path(output, "allele_genotype_summary.tsv"), sep = "\t")
  }
  effect_fields <- c("counted_allele", "difference", "lower", "upper")
  if (is.null(evidence)) {
    hits$counted_allele <- rep(NA_character_, nrow(hits))
    for (field in effect_fields[-1L]) hits[[field]] <- rep(NA_real_, nrow(hits))
  } else hits <- cbind(hits, evidence$estimates[match(hits$marker_id, evidence$estimates$marker_id), effect_fields])
  if (nrow(hits)) {
    support <- stats::aggregate(method ~ marker_id, hits, function(x) paste(methods[methods %in% x], collapse = ", "))
    loci <- unique(hits[c(fields, effect_fields)])
    loci$n_methods <- as.integer(table(hits$marker_id)[loci$marker_id])
    loci$methods <- support$method[match(loci$marker_id, support$marker_id)]
    data.table::fwrite(loci, file.path(output, "chromosome_overview_loci.tsv"), sep = "\t")
  }
  if (!is.null(imputation_accuracy)) files$imputation <- breeding_plot_imputation(imputation_accuracy, output)
  if (!is.null(predictions)) {
    prediction <- breeding_plot_prediction(prediction_input$predictions, observed, prediction_input$ids,
      assessment_type, output, trait, unit, featured, selection_direction)
    files$prediction <- prediction$files
    featured <- prediction$featured
  }
  windows <- data.frame(group_id = character(), chromosome = character(), region_start = numeric(), region_end = numeric(), learner = character())
  if (!is.null(explanation)) {
    regional <- breeding_plot_regional(explanation$importance, hits, !is.null(scans), breeding_plot_genome(maps), output, trait, unit,
                                       explanation$membership)
    windows <- regional$top
    files$regional_importance <- regional$file
    files$local_explanation <- breeding_plot_local(explanation, observed, featured, output, trait, unit)
  }
  if (nrow(hits) || nrow(windows)) {
    overview <- breeding_plot_chromosomes(maps, hits, windows, methods, explanation$receipt$learners, output, trait, unit = unit,
      shap_unit = if (is.null(explanation$membership)) "window" else "group")
    files$chromosome_overview <- overview$files
    data.table::fwrite(overview$pages, file.path(output, "chromosome_overview_pages.tsv"), sep = "\t")
    if (nrow(windows)) {
      overlap <- stats::aggregate(learner ~ group_id + chromosome + region_start + region_end, windows,
        function(x) paste(breeding_plot_learner(explanation$receipt$learners[explanation$receipt$learners %in% x]), collapse = ", "))
      overlap$n_learners <- as.integer(table(windows$group_id)[overlap$group_id])
      overlap$shown_on_chromosome_pages <- overlap$chromosome %in% overview$pages$chromosome
      data.table::fwrite(overlap, file.path(output, "chromosome_overview_ML_windows.tsv"), sep = "\t")
    }
  }
  manifest <- data.frame(figure = rep(names(files), pmax(1L, lengths(files))),
    file = unlist(lapply(files, function(x) if (length(x)) basename(x) else NA_character_), use.names = FALSE))
  manifest$status <- ifelse(is.na(manifest$file), "no_plottable_evidence", "rendered")
  data.table::fwrite(manifest, file.path(output, "figure_manifest.tsv"), sep = "\t", na = "NA")
  list(trait = trait, output = output, figures = unlist(files, use.names = FALSE), tables = tables, manifest = manifest,
       assessment_ids = prediction_input$ids, assessment_id_source = prediction_input$source)
}

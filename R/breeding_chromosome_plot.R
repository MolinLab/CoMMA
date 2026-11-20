breeding_plot_chromosomes <- function(markers, hits, windows, methods, learners, output,
                                      trait, prefix = "08_chromosome_overview", focus = character(), unit = "trait units",
                                      shap_unit = "window", theme = c("dark", "white")) {
  theme <- match.arg(theme)
  empty <- data.frame(chromosome = character(), page = integer())
  if (!nrow(hits) && !nrow(windows)) return(list(files = character(), pages = empty))
  chr <- unique(as.character(markers$chromosome))
  numeric_chr <- suppressWarnings(as.numeric(chr))
  wheat <- grepl("^[1-7][ABD]$", chr)
  chr <- if (all(wheat)) chr[order(substr(chr, 2L, 2L), chr)] else chr[order(numeric_chr, chr)]
  show <- if (nrow(hits)) chr[chr %in% hits$chromosome] else chr[chr %in% windows$chromosome]
  show <- c(intersect(focus, show), setdiff(show, focus))
  pages <- split(show, ceiling(seq_along(show) / 4L))
  if (!length(pages)) stop("No plotted signal matches the chromosome map.")
  keys <- c(methods, learners)
  abbreviations <- c(MLM = "MLM", FarmCPU = "FC", BLINK = "BL", gblup_reml = "GB",
                     elastic_net = "EN", random_forest = "RF", xgboost = "XG", gaussian_rkhs = "RK")
  left <- seq_along(methods) - length(methods) - 1L
  right <- seq_along(learners)
  lanes <- setNames(c(left, right), keys)
  if (anyNA(lanes[c(hits$method, windows$learner)])) stop("Signal learner or method is missing from the tracks.")
  if (anyDuplicated(unique(windows[c("group_id", "chromosome", "region_start", "region_end")])$group_id)) {
    stop("Shared SHAP windows must have identical coordinates.")
  }
  hit_link <- function(rows, id, method, position) {
    if (!nrow(rows)) return(data.frame(chromosome = character(), first = numeric(), last = numeric(), mb = numeric(), n = integer()))
    spans <- lapply(split(rows, rows[[id]]), function(x) {
      pos <- position(x)
      stopifnot(length(unique(x$chromosome)) == 1L, length(unique(pos)) == 1L, !anyDuplicated(x[[method]]))
      positions <- lanes[x[[method]]]
      data.frame(chromosome = x$chromosome[1L], first = min(positions), last = max(positions), mb = pos[1L] / 1e6, n = nrow(x))
    })
    do.call(rbind, spans)
  }
  snps <- hit_link(hits, "marker_id", "method", function(x) x$position_bp)
  groups <- hit_link(windows, "group_id", "learner", function(x) (x$region_start + x$region_end) / 2)
  hits$lane <- unname(lanes[hits$method])
  windows$lane <- unname(lanes[windows$learner])
  hits$direction <- rep("Unknown", nrow(hits))
  known <- is.finite(hits$difference)
  hits$direction[known] <- c("Lower", "Unknown", "Higher")[sign(hits$difference[known]) + 2L]
  hits$open <- known & is.finite(hits$lower) & is.finite(hits$upper) & hits$lower <= 0 & hits$upper >= 0
  if (nrow(hits)) hits <- hits[order(snps$n[match(hits$marker_id, rownames(snps))]), ]
  bg <- "#0D1D2D"; panel <- "#142B3D"; ink <- "#EDF5FA"; muted <- "#AEC5D4"
  gold <- "#FFC069"; mint <- "#60E6C9"; shared_mint <- "#ABE8D8"
  directions <- c(Higher = "#7EC9FF", Lower = "#FF9C91", Unknown = "#D5E0E8")
  rail <- "#355164"; grid_line <- "#294357"; strip <- "#1C384D"
  capsule_edge <- "#65849C"; capsule_fill <- "#B5CAD8"; capsule_highlight <- ink
  if (theme == "white") {
    bg <- "#FFFFFF"; panel <- bg; ink <- "#182C3C"; muted <- "#486171"
    gold <- "#A56315"; mint <- "#00806E"; shared_mint <- "#598277"
    directions <- c(Higher = "#17679D", Lower = "#B74735", Unknown = "#667987")
    rail <- "#B6C9D3"; grid_line <- "#DCE6EB"; strip <- "#EAF2F6"
    capsule_edge <- "#6A8C9F"; capsule_fill <- "#ADC6D3"; capsule_highlight <- "#FFFFFF"
  }
  shapes <- c(Higher = 24, Lower = 25, Unknown = 21)
  spans <- do.call(rbind, lapply(split(markers$position_bp, markers$chromosome), function(x) {
    data.frame(start = min(x) / 1e6, end = max(x) / 1e6)
  }))
  spans$chromosome <- rownames(spans)
  shown_windows <- windows[windows$chromosome %in% show, ]
  max_mb <- max(c(spans$end[spans$chromosome %in% show], shown_windows$region_end / 1e6))
  ticks <- pretty(c(0, max_mb), n = 5)
  max_mb <- max(max_mb, ticks)
  evidence <- if (nrow(hits) && nrow(shown_windows)) "GWAS associations and SHAP regions" else if (nrow(hits))
    "GWAS associations" else "SHAP regions"
  title <- paste(evidence, "by chromosome for", trait)
  stem <- min(c(0, unname(lanes))) - .65
  edge <- max(c(0, unname(lanes))) + .65
  render <- function(chromosomes, page) {
    facet <- function(data) {
      data <- data[data$chromosome %in% chromosomes, , drop = FALSE]
      data$chromosome <- factor(data$chromosome, levels = chromosomes)
      data
    }
    span <- facet(spans)
    gwas <- facet(hits)
    shap <- facet(windows)
    snp_links <- facet(snps[snps$n > 1L, ])
    window_links <- facet(groups[groups$n > 1L, ])
    rails <- merge(span, data.frame(lane = unname(lanes)))
    titles <- setNames(vapply(chromosomes, function(x) {
      if (!nrow(hits)) return(sprintf("%s | %d shared %ss", x, sum(window_links$chromosome == x), shap_unit))
      title <- sprintf("%s | %d GWAS / %d shared", x, sum(gwas$chromosome == x), sum(snp_links$chromosome == x))
      if (length(learners) > 1L) paste0(title, sprintf("\n%d shared ML %ss", sum(window_links$chromosome == x), shap_unit)) else title
    }, character(1L)), chromosomes)
    first <- gwas[order(gwas$p_value, gwas$marker_id, method = "radix"), ]
    first <- first[!duplicated(first$chromosome), ]
    first$estimate <- ifelse(is.finite(first$difference),
                             sprintf("%+.2f %s per %s allele", first$difference, unit, first$counted_allele), "")
    plot <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = rails, ggplot2::aes(x = lane, xend = lane, y = start, yend = end), colour = rail, linewidth = .34) +
      ggplot2::geom_segment(data = span, ggplot2::aes(x = 0, xend = 0, y = start, yend = end), colour = capsule_edge, linewidth = 6.8, lineend = "round") +
      ggplot2::geom_segment(data = span, ggplot2::aes(x = 0, xend = 0, y = start, yend = end), colour = capsule_fill, linewidth = 5.4, lineend = "round") +
      ggplot2::geom_segment(data = span, ggplot2::aes(x = -.07, xend = -.07, y = start, yend = end), colour = capsule_highlight, linewidth = 1, alpha = .55, lineend = "round") +
      ggplot2::geom_segment(data = first, ggplot2::aes(x = lane, xend = -.1, y = position_bp / 1e6, yend = position_bp / 1e6),
                            colour = gold, linewidth = .38, linetype = "dotted") +
      ggplot2::geom_segment(data = snp_links, ggplot2::aes(x = first, xend = last, y = mb, yend = mb), colour = gold, linewidth = .75) +
      ggplot2::geom_segment(data = window_links, ggplot2::aes(x = first, xend = last, y = mb, yend = mb),
                            colour = shared_mint, linetype = "dashed", linewidth = .65) +
      ggplot2::geom_point(data = gwas, ggplot2::aes(x = lane, y = position_bp / 1e6, colour = direction, shape = direction,
                                                   fill = ifelse(open, panel, directions[direction])), size = 2.55, stroke = .50) +
      ggplot2::geom_segment(data = shap, ggplot2::aes(x = lane, xend = lane, y = region_start / 1e6, yend = region_end / 1e6),
                            colour = mint, linewidth = 1.5, lineend = "round") +
      ggplot2::geom_text(data = first, ggplot2::aes(x = (stem + edge) / 2, y = -.081 * max_mb, label = marker_id),
                         colour = ink, size = 3.2, fontface = "bold") +
      ggplot2::geom_text(data = first, ggplot2::aes(x = (stem + edge) / 2, y = -.043 * max_mb, label = estimate, colour = direction), size = 2.9) +
      ggplot2::facet_wrap(~ chromosome, nrow = 1L, labeller = ggplot2::labeller(chromosome = titles)) +
      ggplot2::scale_colour_manual(values = directions, guide = "none") + ggplot2::scale_fill_identity() +
      ggplot2::scale_shape_manual(values = shapes, guide = "none") +
      ggplot2::scale_x_continuous(breaks = unname(lanes), labels = unname(abbreviations[keys]), limits = c(stem, edge)) +
      ggplot2::scale_y_reverse(breaks = ticks, limits = max_mb * c(1.03, -.12)) +
      ggplot2::labs(x = NULL, y = "Position on chromosome (Mb)") + ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = bg, colour = NA),
        panel.background = ggplot2::element_rect(fill = panel, colour = NA), panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(colour = grid_line, linewidth = .26),
        axis.text = ggplot2::element_text(colour = muted, size = 9), axis.title.y = ggplot2::element_text(colour = muted),
        strip.background = ggplot2::element_rect(fill = strip, colour = NA),
        strip.text = ggplot2::element_text(colour = ink, face = "bold", size = 11.5, margin = ggplot2::margin(6, 0, 6, 0)),
        panel.spacing = grid::unit(.19, "in"), plot.margin = ggplot2::margin(3, 24, 4, 23))
    file <- file.path(output, paste0(prefix, if (page > 1L) sprintf("_%02d", page), ".png"))
    grDevices::png(file, width = 13.333, height = 7.5, units = "in", res = 260, type = "cairo", bg = bg)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = bg, col = NA))
    label <- function(text, x, y, color = muted, size = 10, face = "plain", just = "left") {
      grid::grid.text(text, x = x, y = y, just = just, gp = grid::gpar(col = color, fontsize = size, fontface = face))
    }
    label(title, .036, .945, ink, 25, "bold")
    count_label <- function(n, unit) sprintf("%d %s%s", n, unit, if (n == 1L) "" else "s")
    label(sprintf("%s shown  |  %s shared across GWAS methods  |  %s shared across learners",
                  count_label(length(show), "chromosome"), count_label(sum(snps$n > 1L), "SNP"),
                  count_label(sum(groups$n > 1L & groups$chromosome %in% show), shap_unit)), .036, .887, muted, 11.5)
    label(sprintf("%02d / %02d", page, length(pages)), .963, .945, gold, 14, "bold", "right")
    print(plot, newpage = FALSE, vp = grid::viewport(x = .5, y = .482, width = 1, height = .743))
    if (nrow(hits)) {
      for (i in seq_along(c("Higher", "Lower"))) {
        direction <- c("Higher", "Lower")[i]
        x <- c(.042, .157)[i]
        grid::grid.points(grid::unit(x, "npc"), grid::unit(.083, "npc"), pch = shapes[[direction]],
                          size = grid::unit(.105, "in"), gp = grid::gpar(col = directions[[direction]], fill = directions[[direction]]))
        label(paste(direction, "trait"), x + .012, .083, ink)
      }
      grid::grid.segments(.264, .083, .284, .083, gp = grid::gpar(col = gold, lwd = 2))
      label("Same SNP across GWAS methods", .294, .083, ink)
    }
    if (length(learners)) {
      grid::grid.segments(.537, .083, .557, .083, gp = grid::gpar(col = mint, lwd = 3.4))
      label(paste("Top SHAP", shap_unit), .567, .083, ink)
      grid::grid.segments(.721, .083, .747, .083, gp = grid::gpar(col = shared_mint, lwd = 1.7, lty = "dashed"))
      label(paste("Same", shap_unit, "across learners"), .757, .083, ink)
    }
    full <- c(MLM = "MLM", FarmCPU = "FarmCPU", BLINK = "BLINK", gblup_reml = "GBLUP", elastic_net = "Elastic Net",
              random_forest = "Random Forest", xgboost = "XGBoost", gaussian_rkhs = "Gaussian RKHS")
    abbreviations_key <- ifelse(abbreviations[keys] == full[keys], full[keys], paste(abbreviations[keys], full[keys]))
    label(paste(abbreviations_key, collapse = "   |   "), .036, .046, muted, 9)
    note <- if (nrow(hits)) {
      "Open symbols have exploratory 95% intervals spanning zero. Grey indicates a zero or unavailable estimate. Capsules span tested SNPs."
    } else "Capsules span the markers used to fit the learners."
    label(note, .036, .019, muted, 8.7)
    file
  }
  list(files = unlist(Map(render, pages, seq_along(pages)), use.names = FALSE),
       pages = data.frame(chromosome = unlist(pages, use.names = FALSE), page = rep(seq_along(pages), lengths(pages))))
}

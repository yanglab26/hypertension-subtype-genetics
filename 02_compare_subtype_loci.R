#!/usr/bin/env Rscript

# Bidirectional cross-subtype comparison of genome-wide significant loci.
#
# Rule 1: a genome-wide significant SNP in a source locus is in LD (r2 > 0.2)
#         within +/-250 kb with a lead SNP in the target subtype.
# Rule 2: a genome-wide significant SNP in a source locus is within 60 kb of a
#         lead SNP in the target subtype, irrespective of LD.
# Evidence from either direction and either rule defines an overlapping locus
# pair. All lead SNPs are retained; loci with multiple lead SNPs are not reduced
# to a single representative SNP.

suppressPackageStartupMessages(library(data.table))

GWS_P <- 5e-8
LD_R2 <- 0.2
LD_WINDOW_BP <- 250000L
DISTANCE_BP <- 60000L

PHENOTYPE_MAP <- c(
  idh = "IDH", ish = "ISH", sdh = "SDH",
  idhl = "IDHL", ishh = "ISHH",
  IDH = "IDH", ISH = "ISH", SDH = "SDH",
  IDHL = "IDHL", ISHH = "ISHH"
)

normalise_phenotype <- function(x) {
  ans <- unname(PHENOTYPE_MAP[as.character(x)])
  if (anyNA(ans)) stop("Unrecognised phenotype label(s): ", paste(unique(x[is.na(ans)]), collapse = ", "))
  ans
}

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(object_name, " is missing columns: ", paste(missing, collapse = ", "))
}

collapse_unique <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

read_loci <- function(path) {
  x <- fread(path, na.strings = c("NA", "", "."))
  assert_columns(
    x,
    c("phe", "locus_id", "chr", "locus_start", "locus_end", "lead_snps"),
    "Locus table"
  )
  x[, phenotype := normalise_phenotype(phe)]
  x[, `:=`(
    locus_id = as.character(locus_id),
    chr = as.integer(chr),
    locus_start = as.integer(locus_start),
    locus_end = as.integer(locus_end),
    lead_snps = as.character(lead_snps)
  )]
  if (anyDuplicated(x[, .(phenotype, locus_id)])) stop("Duplicate phenotype-locus rows detected.")
  if (x[is.na(chr) | is.na(locus_start) | is.na(locus_end) | locus_start > locus_end, .N]) {
    stop("Invalid locus coordinates detected.")
  }
  x[]
}

expand_lead_snps <- function(loci, gws) {
  leads <- loci[, {
    snps <- trimws(unlist(strsplit(lead_snps, ";", fixed = TRUE)))
    snps <- unique(snps[nzchar(snps) & !is.na(snps)])
    .(lead_snp = snps)
  }, by = .(phenotype, locus_id, chr)]

  annotation <- unique(gws[, .(
    phenotype,
    locus_id,
    lead_snp = study_snp,
    lead_chr = study_chr,
    lead_pos = study_pos,
    lead_p = study_p
  )])
  leads <- merge(leads, annotation, by = c("phenotype", "locus_id", "lead_snp"), all.x = TRUE)
  if (leads[is.na(lead_pos) | is.na(lead_p), .N]) {
    stop("At least one declared lead SNP was not found among its locus-level GWS SNPs.")
  }
  leads[]
}

read_gws <- function(path) {
  x <- fread(path, na.strings = c("NA", "", "."))
  assert_columns(
    x,
    c(
      "study_snp", "phe", "locus_id", "study_chr", "study_pos", "study_p",
      "present_in_EUR2m", "panel_chr", "panel_pos"
    ),
    "Regional GWS-SNP table"
  )
  x[, phenotype := normalise_phenotype(phe)]
  x[, `:=`(
    study_snp = as.character(study_snp),
    locus_id = as.character(locus_id),
    study_chr = as.integer(study_chr),
    study_pos = as.integer(study_pos),
    study_p = as.numeric(study_p),
    present_in_EUR2m = as.logical(present_in_EUR2m),
    panel_chr = as.integer(panel_chr),
    panel_pos = as.integer(panel_pos)
  )]
  if (x[is.na(study_p) | study_p >= GWS_P, .N]) stop("Regional input contains SNPs with P >= 5e-8 or missing P values.")
  unique(x)
}

resolve_manifest_paths <- function(manifest, manifest_path) {
  root <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
  manifest[, path := vapply(path, function(p) {
    p <- as.character(p)
    candidate <- if (grepl("^([A-Za-z]:|/)", p)) p else file.path(root, p)
    normalizePath(candidate, winslash = "/", mustWork = TRUE)
  }, character(1L))]
  manifest
}

read_plink_ld <- function(path, phenotype) {
  x <- fread(path, showProgress = FALSE)
  setnames(x, toupper(trimws(names(x))))
  assert_columns(x, c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"), basename(path))
  x <- x[, .(
    phenotype = phenotype,
    query_chr = as.integer(CHR_A),
    query_pos = as.integer(BP_A),
    query_snp = as.character(SNP_A),
    proxy_chr = as.integer(CHR_B),
    proxy_pos = as.integer(BP_B),
    proxy_snp = as.character(SNP_B),
    r2 = as.numeric(R2)
  )]
  if (x[is.na(r2) | r2 < 0 | r2 > 1, .N]) stop("Invalid r2 values in ", path)
  x[, distance_query_proxy := abs(proxy_pos - query_pos)]
  unique(x[distance_query_proxy <= LD_WINDOW_BP & r2 > LD_R2])
}

read_ld_manifest <- function(path, gws) {
  manifest <- fread(path)
  assert_columns(manifest, c("phenotype", "path"), "LD manifest")
  manifest[, phenotype := normalise_phenotype(phenotype)]
  manifest <- resolve_manifest_paths(manifest, path)
  ld <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
    read_plink_ld(manifest$path[i], manifest$phenotype[i])
  }), use.names = TRUE)

  # PLINK may omit self-pairs; add them only for reference-panel-evaluable SNPs.
  self <- unique(gws[present_in_EUR2m == TRUE, .(
    phenotype,
    query_chr = panel_chr,
    query_pos = panel_pos,
    query_snp = study_snp,
    proxy_chr = panel_chr,
    proxy_pos = panel_pos,
    proxy_snp = study_snp,
    r2 = 1,
    distance_query_proxy = 0L
  )])
  unique(rbind(ld, self, fill = TRUE))
}

directional_compare <- function(source, target, gws, leads, ld) {
  source_gws <- gws[phenotype == source, .(
    source_locus = locus_id,
    source_gws_snp = study_snp,
    source_chr = study_chr,
    source_pos = study_pos,
    source_p = study_p
  )]
  target_leads <- leads[phenotype == target, .(
    target_locus = locus_id,
    target_lead_snp = lead_snp,
    target_chr = lead_chr,
    target_pos = lead_pos,
    target_p = lead_p
  )]

  source_ld <- ld[phenotype == source, .(
    source_gws_snp = query_snp,
    ld_query_chr = query_chr,
    ld_query_pos = query_pos,
    target_lead_snp = proxy_snp,
    ld_proxy_chr = proxy_chr,
    ld_proxy_pos = proxy_pos,
    r2,
    distance_query_proxy
  )]
  ld_hits <- merge(source_gws, source_ld, by = "source_gws_snp", allow.cartesian = TRUE)
  ld_hits <- ld_hits[source_chr == ld_query_chr & source_pos == ld_query_pos & r2 > LD_R2 & distance_query_proxy <= LD_WINDOW_BP]
  ld_hits <- merge(ld_hits, target_leads, by = "target_lead_snp", allow.cartesian = TRUE)
  ld_hits <- ld_hits[ld_proxy_chr == target_chr & ld_proxy_pos == target_pos]
  if (nrow(ld_hits)) {
    ld_hits <- ld_hits[, .(
      source_phenotype = source, source_locus, source_gws_snp,
      source_chr, source_pos, source_p,
      target_phenotype = target, target_locus, target_lead_snp,
      target_chr, target_pos, target_p,
      evidence_rule = "LD_r2_gt_0.2_within_250kb",
      r2, distance_bp = abs(source_pos - target_pos)
    )]
  } else ld_hits <- data.table()

  distance_hits <- merge(source_gws, target_leads, by.x = "source_chr", by.y = "target_chr", allow.cartesian = TRUE)
  setnames(distance_hits, "source_chr", "shared_chr")
  distance_hits[, distance_bp := abs(source_pos - target_pos)]
  distance_hits <- distance_hits[distance_bp <= DISTANCE_BP]
  if (nrow(distance_hits)) {
    distance_hits <- distance_hits[, .(
      source_phenotype = source, source_locus, source_gws_snp,
      source_chr = shared_chr, source_pos, source_p,
      target_phenotype = target, target_locus, target_lead_snp,
      target_chr = shared_chr, target_pos, target_p,
      evidence_rule = "distance_le_60kb", r2 = NA_real_, distance_bp
    )]
  } else distance_hits <- data.table()

  out <- unique(rbind(ld_hits, distance_hits, fill = TRUE))
  if (nrow(out)) out[, direction := paste0(source, "_to_", target)]
  out
}

compare_pair <- function(a, b, gws, leads, loci, ld) {
  evidence <- unique(rbind(
    directional_compare(a, b, gws, leads, ld),
    directional_compare(b, a, gws, leads, ld),
    fill = TRUE
  ))
  if (nrow(evidence)) {
    evidence[, `:=`(
      phenotype_A = a,
      phenotype_B = b,
      locus_A = fifelse(source_phenotype == a, source_locus, target_locus),
      locus_B = fifelse(source_phenotype == b, source_locus, target_locus)
    )]
    pairs <- evidence[, .(
      n_supporting_rows = .N,
      directions = collapse_unique(direction),
      evidence_rules = collapse_unique(evidence_rule),
      supported_by_LD = any(evidence_rule == "LD_r2_gt_0.2_within_250kb"),
      supported_by_60kb = any(evidence_rule == "distance_le_60kb"),
      maximum_r2 = safe_max(r2),
      minimum_distance_bp = safe_min(distance_bp),
      source_GWS_SNPs = collapse_unique(paste(source_phenotype, source_gws_snp, sep = ":")),
      target_lead_SNPs = collapse_unique(paste(target_phenotype, target_lead_snp, sep = ":"))
    ), by = .(phenotype_A, locus_A, phenotype_B, locus_B)]
  } else {
    pairs <- data.table(
      phenotype_A = character(), locus_A = character(), phenotype_B = character(), locus_B = character(),
      n_supporting_rows = integer(), directions = character(), evidence_rules = character(),
      supported_by_LD = logical(), supported_by_60kb = logical(), maximum_r2 = numeric(),
      minimum_distance_bp = numeric(), source_GWS_SNPs = character(), target_lead_SNPs = character()
    )
  }

  status <- rbind(
    loci[phenotype == a, .(phenotype, locus_id, chr, locus_start, locus_end,
      comparator = b, overlap = locus_id %in% pairs$locus_A)],
    loci[phenotype == b, .(phenotype, locus_id, chr, locus_start, locus_end,
      comparator = a, overlap = locus_id %in% pairs$locus_B)]
  )
  list(evidence = evidence, pairs = pairs, status = status)
}

build_strict_triplets <- function(pair_results) {
  idh_ish <- unique(pair_results[["IDH_vs_ISH"]]$pairs[, .(idh_locus = locus_A, ish_locus = locus_B)])
  idh_sdh <- unique(pair_results[["IDH_vs_SDH"]]$pairs[, .(idh_locus = locus_A, sdh_locus = locus_B)])
  ish_sdh <- unique(pair_results[["ISH_vs_SDH"]]$pairs[, .(ish_locus = locus_A, sdh_locus = locus_B)])
  candidates <- merge(idh_ish, idh_sdh, by = "idh_locus", allow.cartesian = TRUE)
  unique(merge(candidates, ish_sdh, by = c("ish_locus", "sdh_locus")))
}

run_cross_subtype_comparison <- function(locus_file, gws_file, ld_manifest_file, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  loci <- read_loci(locus_file)
  gws <- read_gws(gws_file)
  leads <- expand_lead_snps(loci, gws)
  ld <- read_ld_manifest(ld_manifest_file, gws)

  available <- unique(loci$phenotype)
  planned <- list(c("IDH", "ISH"), c("IDH", "SDH"), c("ISH", "SDH"), c("IDH", "IDHL"), c("ISH", "ISHH"))
  planned <- Filter(function(z) all(z %in% available), planned)

  results <- list()
  for (pair in planned) {
    key <- paste(pair, collapse = "_vs_")
    results[[key]] <- compare_pair(pair[1L], pair[2L], gws, leads, loci, ld)
  }

  all_evidence <- rbindlist(lapply(names(results), function(k) {
    x <- copy(results[[k]]$evidence); if (nrow(x)) x[, comparison := k]; x
  }), fill = TRUE)
  all_pairs <- rbindlist(lapply(names(results), function(k) {
    x <- copy(results[[k]]$pairs); if (nrow(x)) x[, comparison := k]; x
  }), fill = TRUE)
  all_status <- rbindlist(lapply(names(results), function(k) {
    x <- copy(results[[k]]$status); x[, comparison := k]; x
  }), fill = TRUE)

  fwrite(all_evidence, file.path(output_dir, "01_bidirectional_SNP_level_evidence.tsv"), sep = "\t", na = "NA")
  fwrite(all_pairs, file.path(output_dir, "02_cross_subtype_locus_pairs.tsv"), sep = "\t", na = "NA")
  fwrite(all_status, file.path(output_dir, "03_pairwise_locus_status.tsv"), sep = "\t", na = "NA")
  fwrite(leads, file.path(output_dir, "04_expanded_lead_SNP_manifest.tsv"), sep = "\t", na = "NA")

  if (all(c("IDH_vs_ISH", "IDH_vs_SDH", "ISH_vs_SDH") %in% names(results))) {
    triplets <- build_strict_triplets(results)
    main_edges <- rbindlist(list(
      results[["IDH_vs_ISH"]]$pairs[, .(phenotype_A, locus_A, phenotype_B, locus_B)],
      results[["IDH_vs_SDH"]]$pairs[, .(phenotype_A, locus_A, phenotype_B, locus_B)],
      results[["ISH_vs_SDH"]]$pairs[, .(phenotype_A, locus_A, phenotype_B, locus_B)]
    ))
    fwrite(triplets, file.path(output_dir, "05_strict_three_subtype_locus_triplets.tsv"), sep = "\t", na = "NA")
    fwrite(main_edges, file.path(output_dir, "06_main_subtype_pairwise_edges.tsv"), sep = "\t", na = "NA")
  }

  parameters <- data.table(
    parameter = c("genome_build", "GWS_P", "LD_reference", "LD_r2_rule", "LD_window_bp", "distance_rule_bp"),
    value = c("GRCh37", format(GWS_P, scientific = TRUE), "1000 Genomes Phase 3 European", "r2 > 0.2", LD_WINDOW_BP, DISTANCE_BP)
  )
  fwrite(parameters, file.path(output_dir, "00_analysis_parameters.tsv"), sep = "\t")
  invisible(results)
}

run_locus_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) {
    stop("Usage: Rscript 02_compare_subtype_loci.R <loci.tsv> <regional_gws_snps.tsv> <ld_manifest.tsv> <output_dir>")
  }
  run_cross_subtype_comparison(args[1L], args[2L], args[3L], args[4L])
}

if (sys.nframe() == 0L) run_locus_cli()

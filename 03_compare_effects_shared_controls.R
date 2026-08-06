#!/usr/bin/env Rscript

# Shared-control-adjusted comparison of subtype GWAS effect estimates.
#
# The Lin-Sullivan sampling correlation is calculated from the case/control
# sample composition and is used in a covariance-adjusted Delta/Wald test of
# beta_1 - beta_2 = 0. The primary multiplicity family contains every unique
# lead SNP crossed with all three pairwise subtype comparisons.
#
# Reference: Lin DY, Sullivan PF. Am J Hum Genet. 2009;85:862-872.
# PMID: 20004761.

suppressPackageStartupMessages(library(data.table))

ALPHA <- 0.05
GWS_P <- 5e-8
PRIMARY_PHENOTYPES <- c("IDH", "ISH", "SDH")

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(object_name, " is missing columns: ", paste(missing, collapse = ", "))
}

collapse_unique <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

resolve_paths <- function(manifest, manifest_path) {
  root <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
  manifest[, path := vapply(path, function(p) {
    p <- as.character(p)
    candidate <- if (grepl("^([A-Za-z]:|/)", p)) p else file.path(root, p)
    normalizePath(candidate, winslash = "/", mustWork = TRUE)
  }, character(1L))]
  manifest
}

read_primary_loci <- function(path) {
  x <- fread(path, na.strings = c("NA", "", "."))
  assert_columns(
    x,
    c("phe", "genomic_locus", "locus_id", "chr", "locus_start", "locus_end", "top_snp", "top_pos", "top_p", "n_lead_snps", "lead_snps"),
    "Locus table"
  )
  x[, phenotype := toupper(as.character(phe))]
  x <- x[phenotype %chin% PRIMARY_PHENOTYPES]
  x[, `:=`(
    locus_id = as.character(locus_id),
    chr = as.integer(chr),
    locus_start = as.integer(locus_start),
    locus_end = as.integer(locus_end),
    top_snp = as.character(top_snp),
    top_pos = as.integer(top_pos),
    top_p = as.numeric(top_p),
    n_lead_snps = as.integer(n_lead_snps),
    lead_snps = as.character(lead_snps)
  )]
  if (anyDuplicated(x[, .(phenotype, locus_id)])) stop("Duplicate primary phenotype-locus records detected.")
  x[]
}

expand_lead_manifest <- function(loci) {
  manifest <- loci[, {
    snps <- trimws(unlist(strsplit(lead_snps, ";", fixed = TRUE)))
    snps <- snps[nzchar(snps) & !is.na(snps)]
    .(SNP = snps)
  }, by = .(
    phenotype, genomic_locus, locus_id, locus_chr = chr, locus_start, locus_end,
    top_snp, top_pos, top_p, declared_n_lead_snps = n_lead_snps
  )]
  manifest[, is_top_snp := SNP == top_snp]

  qc <- manifest[, .(
    expanded_n_lead_snps = .N,
    declared_n_lead_snps = unique(declared_n_lead_snps)
  ), by = .(phenotype, locus_id)]
  if (qc[lengths(declared_n_lead_snps) != 1L | expanded_n_lead_snps != declared_n_lead_snps, .N]) {
    stop("Expanded lead-SNP counts do not agree with n_lead_snps.")
  }
  unique(manifest)
}

read_gwas_design <- function(path) {
  x <- fread(path)
  assert_columns(x, c("phenotype", "path", "n_case", "n_control"), "GWAS manifest")
  x[, phenotype := toupper(as.character(phenotype))]
  if (!setequal(x$phenotype, PRIMARY_PHENOTYPES) || nrow(x) != 3L) {
    stop("GWAS manifest must contain exactly one row each for IDH, ISH, and SDH.")
  }
  x[, `:=`(n_case = as.numeric(n_case), n_control = as.numeric(n_control))]
  if (x[!is.finite(n_case) | !is.finite(n_control) | n_case <= 0 | n_control <= 0, .N]) stop("Invalid sample sizes in GWAS manifest.")
  resolve_paths(x, path)
}

read_overlap_design <- function(path, gwas_design) {
  x <- fread(path)
  assert_columns(x, c("phenotype_1", "phenotype_2", "n_overlap_cases", "n_overlap_controls"), "Overlap manifest")
  x[, `:=`(
    phenotype_1 = toupper(as.character(phenotype_1)),
    phenotype_2 = toupper(as.character(phenotype_2)),
    n_overlap_cases = as.numeric(n_overlap_cases),
    n_overlap_controls = as.numeric(n_overlap_controls)
  )]
  expected <- c("IDH_vs_ISH", "IDH_vs_SDH", "ISH_vs_SDH")
  x[, comparison := paste(phenotype_1, phenotype_2, sep = "_vs_")]
  if (!setequal(x$comparison, expected) || nrow(x) != 3L) {
    stop("Overlap manifest must contain IDH-vs-ISH, IDH-vs-SDH, and ISH-vs-SDH once each.")
  }
  x <- merge(x, gwas_design[, .(phenotype_1 = phenotype, n_case_1 = n_case, n_control_1 = n_control)], by = "phenotype_1")
  x <- merge(x, gwas_design[, .(phenotype_2 = phenotype, n_case_2 = n_case, n_control_2 = n_control)], by = "phenotype_2")
  x[]
}

read_one_gwas <- function(path, phenotype, target_snps) {
  required <- c("SNP", "chr", "pos", "A1", "A2", "beta", "se", "OR", "eaf", "pval", "maf")
  header <- names(fread(path, nrows = 0L))
  missing <- setdiff(required, header)
  if (length(missing)) stop(phenotype, " GWAS is missing columns: ", paste(missing, collapse = ", "))
  x <- fread(path, select = required, showProgress = FALSE)
  x <- x[SNP %chin% target_snps]
  x[, `:=`(
    phenotype = phenotype,
    SNP = as.character(SNP), chr = as.character(chr), pos = as.integer(pos),
    A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)),
    beta = as.numeric(beta), se = as.numeric(se), OR = as.numeric(OR),
    eaf = as.numeric(eaf), pval = as.numeric(pval), maf = as.numeric(maf)
  )]
  if (x[, anyDuplicated(SNP)]) stop(phenotype, " GWAS contains duplicate target rsIDs.")
  x[]
}

complement_allele <- function(x) chartr("ACGT", "TGCA", x)

is_snv_allele <- function(x) !is.na(x) & nchar(x) == 1L & x %chin% c("A", "C", "G", "T")

is_palindromic <- function(a1, a2) paste0(a1, a2) %chin% c("AT", "TA", "CG", "GC")

harmonize_effects <- function(gwas_long, target_snps) {
  expected <- CJ(SNP = target_snps, phenotype = PRIMARY_PHENOTYPES, unique = TRUE)
  found <- unique(gwas_long[, .(SNP, phenotype, found = is.finite(beta) & is.finite(se) & se > 0)])
  extraction_qc <- merge(expected, found, by = c("SNP", "phenotype"), all.x = TRUE)
  extraction_qc[is.na(found), found := FALSE]

  coordinate_qc <- gwas_long[, .(
    n_distinct_chr_pos = uniqueN(paste(chr, pos, sep = ":")),
    chr_pos_values = collapse_unique(paste(chr, pos, sep = ":"))
  ), by = SNP]

  gwas_long[, phenotype_order := match(phenotype, PRIMARY_PHENOTYPES)]
  reference <- gwas_long[order(SNP, phenotype_order), .(
    ref_A1 = A1[1L], ref_A2 = A2[1L], reference_phenotype = phenotype[1L]
  ), by = SNP]
  h <- merge(gwas_long, reference, by = "SNP", all.x = TRUE)
  h <- merge(h, coordinate_qc, by = "SNP", all.x = TRUE)
  h[, simple_snv := is_snv_allele(A1) & is_snv_allele(A2) & is_snv_allele(ref_A1) & is_snv_allele(ref_A2)]
  h[, palindromic := is_palindromic(A1, A2)]
  h[, allele_status := fcase(
    A1 == ref_A1 & A2 == ref_A2, "direct",
    !palindromic & A1 == ref_A2 & A2 == ref_A1, "swapped",
    simple_snv & !palindromic & complement_allele(A1) == ref_A1 & complement_allele(A2) == ref_A2, "strand_direct",
    simple_snv & !palindromic & complement_allele(A1) == ref_A2 & complement_allele(A2) == ref_A1, "strand_swapped",
    palindromic & A1 == ref_A2 & A2 == ref_A1, "palindromic_swapped_ambiguous",
    default = "mismatch"
  )]
  h[, valid_harmonization := allele_status %chin% c("direct", "swapped", "strand_direct", "strand_swapped") &
      is.finite(beta) & is.finite(se) & se > 0 & n_distinct_chr_pos == 1L]
  h[, flip := allele_status %chin% c("swapped", "strand_swapped")]
  h[, `:=`(
    effect_allele = ref_A1,
    other_allele = ref_A2,
    beta_harmonized = fifelse(flip, -beta, beta),
    se_harmonized = se,
    eaf_harmonized = fifelse(flip, 1 - eaf, eaf)
  )]

  qc <- h[, .(SNP, phenotype, chr, pos, A1, A2, ref_A1, ref_A2, allele_status, valid_harmonization, n_distinct_chr_pos, chr_pos_values)]
  list(
    data = h[valid_harmonization == TRUE, .(
      SNP, phenotype, chr, pos, effect_allele, other_allele,
      beta = beta_harmonized, se = se_harmonized, OR = exp(beta_harmonized),
      eaf = eaf_harmonized, pval, maf, allele_status, original_A1 = A1, original_A2 = A2
    )],
    extraction_qc = extraction_qc,
    harmonization_qc = qc,
    complete = extraction_qc[found == FALSE, .N] == 0L && qc[valid_harmonization == FALSE, .N] == 0L
  )
}

lin_sullivan_rho <- function(n_case_1, n_control_1, n_case_2, n_control_2,
                             n_overlap_cases = 0, n_overlap_controls = 0) {
  values <- as.numeric(c(n_case_1, n_control_1, n_case_2, n_control_2, n_overlap_cases, n_overlap_controls))
  n_case_1 <- values[1L]; n_control_1 <- values[2L]
  n_case_2 <- values[3L]; n_control_2 <- values[4L]
  n_overlap_cases <- values[5L]; n_overlap_controls <- values[6L]
  numerator <- n_overlap_controls * sqrt((n_case_1 * n_case_2) / (n_control_1 * n_control_2)) +
    n_overlap_cases * sqrt((n_control_1 * n_control_2) / (n_case_1 * n_case_2))
  numerator / sqrt((n_case_1 + n_control_1) * (n_case_2 + n_control_2))
}

two_sided_normal_p <- function(z) {
  log_p <- log(2) + pnorm(-abs(z), log.p = TRUE)
  pmax(exp(log_p), .Machine$double.xmin)
}

compare_one_pair <- function(p1, p2, rho, effects) {
  x1 <- effects[phenotype == p1, .(
    SNP, chr, pos, effect_allele, other_allele,
    beta_1 = beta, se_1 = se, OR_1 = OR, eaf_1 = eaf, pval_1 = pval
  )]
  x2 <- effects[phenotype == p2, .(
    SNP, beta_2 = beta, se_2 = se, OR_2 = OR, eaf_2 = eaf, pval_2 = pval
  )]
  ans <- merge(x1, x2, by = "SNP")
  ans[, `:=`(
    comparison = paste(p1, p2, sep = "_vs_"), phenotype_1 = p1, phenotype_2 = p2,
    sampling_correlation_rho = rho,
    covariance_beta_1_beta_2 = rho * se_1 * se_2
  )]
  ans[, variance_difference := se_1^2 + se_2^2 - 2 * covariance_beta_1_beta_2]
  if (ans[!is.finite(variance_difference) | variance_difference <= 0, .N]) stop("Invalid adjusted variance for ", p1, " vs ", p2)
  ans[, `:=`(
    beta_difference = beta_1 - beta_2,
    se_difference_adjusted = sqrt(variance_difference)
  )]
  ans[, z_difference_adjusted := beta_difference / se_difference_adjusted]
  ans[, `:=`(
    p_difference_adjusted = two_sided_normal_p(z_difference_adjusted),
    ratio_of_ORs = exp(beta_difference),
    effect_direction = fifelse(sign(beta_1) == sign(beta_2), "same_direction", "opposite_direction"),
    gws_1 = pval_1 < GWS_P,
    gws_2 = pval_2 < GWS_P
  )]
  ans[]
}

build_region_map <- function(loci, edges) {
  assert_columns(edges, c("locus_A", "locus_B"), "Cross-subtype edge table")
  ids <- unique(loci$locus_id)
  adjacency <- setNames(vector("list", length(ids)), ids)
  for (i in seq_len(nrow(edges))) {
    a <- as.character(edges$locus_A[i]); b <- as.character(edges$locus_B[i])
    if (a %chin% ids && b %chin% ids) {
      adjacency[[a]] <- unique(c(adjacency[[a]], b)); adjacency[[b]] <- unique(c(adjacency[[b]], a))
    }
  }
  visited <- setNames(rep(FALSE, length(ids)), ids)
  components <- list(); k <- 0L
  for (start in ids) {
    if (visited[[start]]) next
    k <- k + 1L; queue <- start; members <- character()
    while (length(queue)) {
      current <- queue[1L]; queue <- queue[-1L]
      if (visited[[current]]) next
      visited[[current]] <- TRUE; members <- c(members, current)
      queue <- c(queue, adjacency[[current]][!visited[adjacency[[current]]]])
    }
    components[[k]] <- data.table(locus_id = members, component = k)
  }
  map <- merge(loci, rbindlist(components), by = "locus_id")
  region <- map[, .(
    region_chr = chr[which.min(locus_start)],
    region_start = min(locus_start), region_end = max(locus_end),
    member_phenotypes = collapse_unique(phenotype), member_loci = collapse_unique(locus_id)
  ), by = component]
  region[, chr_order := suppressWarnings(as.integer(gsub("^chr", "", as.character(region_chr), ignore.case = TRUE)))]
  setorder(region, chr_order, region_start)
  region[, unified_region_id := paste0("Region_", seq_len(.N))]
  region[, chr_order := NULL]
  merge(map, region, by = "component")
}

run_shared_control_comparison <- function(locus_file, gwas_manifest_file, overlap_manifest_file,
                                          edge_file = NA_character_, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  loci <- read_primary_loci(locus_file)
  lead_manifest <- expand_lead_manifest(loci)
  unique_leads <- lead_manifest[, .(
    source_phenotypes = collapse_unique(phenotype),
    source_loci = collapse_unique(locus_id),
    n_source_records = .N,
    is_top_snp_any = any(is_top_snp)
  ), by = SNP]

  design <- read_gwas_design(gwas_manifest_file)
  overlap <- read_overlap_design(overlap_manifest_file, design)
  gwas_long <- rbindlist(lapply(PRIMARY_PHENOTYPES, function(ph) {
    read_one_gwas(design[phenotype == ph, path], ph, unique_leads$SNP)
  }), use.names = TRUE)

  h <- harmonize_effects(gwas_long, unique_leads$SNP)
  fwrite(h$extraction_qc, file.path(output_dir, "03a_extraction_QC.tsv"), sep = "\t", na = "NA")
  fwrite(h$harmonization_qc, file.path(output_dir, "03b_harmonization_QC.tsv"), sep = "\t", na = "NA")
  if (!h$complete) {
    stop("Incomplete extraction or invalid/ambiguous allele harmonization; inspect 03a and 03b QC outputs.")
  }
  effects <- h$data

  overlap[, sampling_correlation_rho := mapply(
    lin_sullivan_rho, n_case_1, n_control_1, n_case_2, n_control_2,
    n_overlap_cases, n_overlap_controls
  )]
  if (overlap[!is.finite(sampling_correlation_rho) | abs(sampling_correlation_rho) >= 1, .N]) {
    stop("Invalid Lin-Sullivan sampling correlation.")
  }
  rho_matrix <- diag(3L); dimnames(rho_matrix) <- list(PRIMARY_PHENOTYPES, PRIMARY_PHENOTYPES)
  for (i in seq_len(nrow(overlap))) {
    rho_matrix[overlap$phenotype_1[i], overlap$phenotype_2[i]] <- overlap$sampling_correlation_rho[i]
    rho_matrix[overlap$phenotype_2[i], overlap$phenotype_1[i]] <- overlap$sampling_correlation_rho[i]
  }
  if (min(eigen(rho_matrix, symmetric = TRUE, only.values = TRUE)$values) <= 0) stop("Sampling-correlation matrix is not positive definite.")

  comparisons <- rbindlist(lapply(seq_len(nrow(overlap)), function(i) {
    compare_one_pair(overlap$phenotype_1[i], overlap$phenotype_2[i], overlap$sampling_correlation_rho[i], effects)
  }))
  expected_tests <- nrow(unique_leads) * 3L
  if (nrow(comparisons) != expected_tests) stop("Unexpected number of formal pairwise tests.")
  comparisons <- merge(comparisons, unique_leads, by = "SNP", all.x = TRUE)
  threshold <- ALPHA / nrow(comparisons)
  comparisons[, `:=`(
    p_bonferroni_global = p.adjust(p_difference_adjusted, method = "bonferroni"),
    significant_bonferroni_global = p_difference_adjusted < threshold
  )]
  setorder(comparisons, comparison, p_difference_adjusted, SNP)
  significant <- comparisons[significant_bonferroni_global == TRUE]
  significant[, `:=`(
    beta_1_lower_95 = beta_1 - qnorm(0.975) * se_1,
    beta_1_upper_95 = beta_1 + qnorm(0.975) * se_1,
    beta_2_lower_95 = beta_2 - qnorm(0.975) * se_2,
    beta_2_upper_95 = beta_2 + qnorm(0.975) * se_2
  )]

  source_results <- merge(lead_manifest, comparisons, by = "SNP", allow.cartesian = TRUE)
  fwrite(lead_manifest, file.path(output_dir, "01_lead_SNP_source_manifest.tsv"), sep = "\t")
  fwrite(unique_leads, file.path(output_dir, "02_unique_lead_SNPs.tsv"), sep = "\t")
  fwrite(effects, file.path(output_dir, "04_harmonized_effects.tsv"), sep = "\t")
  fwrite(overlap, file.path(output_dir, "05_shared_control_design.tsv"), sep = "\t")
  fwrite(comparisons, file.path(output_dir, "06_all_effect_comparisons.tsv"), sep = "\t")
  fwrite(significant, file.path(output_dir, "07_global_Bonferroni_significant.tsv"), sep = "\t")
  fwrite(significant[, .(
    SNP, chr, pos, effect_allele, other_allele, comparison, phenotype_1, phenotype_2,
    beta_1, beta_1_lower_95, beta_1_upper_95, beta_2, beta_2_lower_95, beta_2_upper_95,
    beta_difference, se_difference_adjusted, p_difference_adjusted
  )], file.path(output_dir, "08_significant_forest_plot_data.tsv"), sep = "\t")
  fwrite(source_results, file.path(output_dir, "09_all_comparisons_with_source_loci.tsv"), sep = "\t")

  n_regions <- NA_integer_
  if (!is.na(edge_file) && nzchar(edge_file) && toupper(edge_file) != "NONE") {
    edges <- fread(edge_file)
    region_map <- build_region_map(loci, edges)
    mapped <- merge(
      significant,
      unique(lead_manifest[, .(SNP, source_locus = locus_id)]),
      by = "SNP", allow.cartesian = TRUE
    )
    mapped <- merge(mapped, region_map[, .(source_locus = locus_id, unified_region_id, region_chr, region_start, region_end)], by = "source_locus")
    region_summary <- unique(mapped, by = c("comparison", "SNP", "unified_region_id"))[, .(
      significant_SNPs = collapse_unique(SNP),
      n_significant_SNPs = uniqueN(SNP),
      minimum_P_difference = min(p_difference_adjusted)
    ), by = .(comparison, unified_region_id, region_chr, region_start, region_end)]
    n_regions <- uniqueN(region_summary$unified_region_id)
    fwrite(region_map, file.path(output_dir, "10_unified_region_map.tsv"), sep = "\t")
    fwrite(region_summary, file.path(output_dir, "11_significant_unified_regions.tsv"), sep = "\t")
  }

  summary <- data.table(
    metric = c("source_loci", "lead_SNP_source_records", "unique_lead_SNPs", "formal_tests", "global_Bonferroni_threshold", "significant_tests", "significant_unique_SNPs", "significant_unified_regions"),
    value = c(nrow(loci), nrow(lead_manifest), nrow(unique_leads), nrow(comparisons), threshold, nrow(significant), uniqueN(significant$SNP), n_regions)
  )
  fwrite(summary, file.path(output_dir, "00_analysis_summary.tsv"), sep = "\t", na = "NA")
  invisible(list(comparisons = comparisons, significant = significant, summary = summary))
}

run_effect_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 5L) {
    stop("Usage: Rscript 03_compare_effects_shared_controls.R <loci.tsv> <gwas_manifest.tsv> <overlap_manifest.tsv> <cross_subtype_edges.tsv|NONE> <output_dir>")
  }
  run_shared_control_comparison(args[1L], args[2L], args[3L], args[4L], args[5L])
}

if (sys.nframe() == 0L) run_effect_cli()

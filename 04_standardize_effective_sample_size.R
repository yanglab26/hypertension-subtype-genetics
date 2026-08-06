#!/usr/bin/env Rscript

# Effective-sample-size standardization of case-control GWAS summary statistics.
#
# N_eff = 4 / (1 / N_case + 1 / N_control)
# SE_standardized = SE_original * sqrt(N_eff_original / N_eff_target)
#
# Effect estimates are unchanged. This is a summary-statistic sensitivity
# analysis and does not resample participants or rerun the GWAS.

suppressPackageStartupMessages(library(data.table))

GWS_P <- 5e-8
PRIMARY_PHENOTYPES <- c("IDH", "ISH", "SDH")

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(object_name, " is missing columns: ", paste(missing, collapse = ", "))
}

two_sided_normal_p <- function(z) {
  log_p <- log(2) + pnorm(-abs(z), log.p = TRUE)
  pmax(exp(log_p), .Machine$double.xmin)
}

read_sample_design <- function(path) {
  x <- fread(path)
  assert_columns(x, c("phenotype", "n_case", "n_control"), "Sample-size table")
  x[, phenotype := toupper(as.character(phenotype))]
  if (!setequal(x$phenotype, PRIMARY_PHENOTYPES) || nrow(x) != 3L) {
    stop("Sample-size table must contain exactly one row each for IDH, ISH, and SDH.")
  }
  x[, `:=`(n_case = as.numeric(n_case), n_control = as.numeric(n_control))]
  if (x[!is.finite(n_case) | !is.finite(n_control) | n_case <= 0 | n_control <= 0, .N]) {
    stop("Invalid case/control sample sizes.")
  }
  x[, n_eff := 4 / (1 / n_case + 1 / n_control)]
  target <- min(x$n_eff)
  x[, `:=`(
    target_n_eff = target,
    information_fraction_retained = target / n_eff,
    se_scale_factor = sqrt(n_eff / target)
  )]
  x[]
}

standardize_effects <- function(effects, sample_design) {
  assert_columns(
    effects,
    c("SNP", "phenotype", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval"),
    "Harmonized effect table"
  )
  effects[, phenotype := toupper(as.character(phenotype))]
  effects <- effects[phenotype %chin% PRIMARY_PHENOTYPES]
  effects[, `:=`(beta = as.numeric(beta), se = as.numeric(se), pval = as.numeric(pval))]
  if (effects[!is.finite(beta) | !is.finite(se) | se <= 0, .N]) stop("Invalid beta or SE values.")
  if (anyDuplicated(effects[, .(SNP, phenotype)])) stop("Duplicate SNP-phenotype rows detected.")

  x <- merge(effects, sample_design, by = "phenotype", all.x = TRUE)
  if (x[is.na(n_eff), .N]) stop("A phenotype is missing from the sample-size table.")
  x[, `:=`(
    z_original_recalculated = beta / se,
    p_original_recalculated = two_sided_normal_p(beta / se),
    se_standardized = se * se_scale_factor
  )]
  x[, `:=`(
    z_standardized = beta / se_standardized,
    p_standardized = two_sided_normal_p(beta / se_standardized)
  )]
  x[, `:=`(
    gws_original = pval < GWS_P,
    gws_standardized = p_standardized < GWS_P
  )]
  x[, significance_transition := fcase(
    gws_original & gws_standardized, "GWS_retained",
    gws_original & !gws_standardized, "GWS_not_retained",
    !gws_original & gws_standardized, "became_GWS_after_standardization",
    default = "not_GWS"
  )]
  x[]
}

summarize_source_loci <- function(standardized, lead_manifest) {
  assert_columns(
    lead_manifest,
    c("phenotype", "genomic_locus", "locus_id", "locus_chr", "locus_start", "locus_end", "SNP", "is_top_snp"),
    "Lead-SNP source manifest"
  )
  lead_manifest[, phenotype := toupper(as.character(phenotype))]
  source <- merge(
    lead_manifest,
    standardized[, .(
      phenotype, SNP, beta, se_original = se, p_original = pval,
      se_standardized, p_standardized, gws_original, gws_standardized,
      n_eff, target_n_eff, se_scale_factor, information_fraction_retained
    )],
    by = c("phenotype", "SNP"),
    all.x = TRUE
  )
  if (source[is.na(p_standardized), .N]) stop("At least one source lead SNP lacks a standardized result in its discovery phenotype.")

  locus <- source[, {
    best <- which.min(p_standardized)
    .(
      n_predefined_lead_SNPs = .N,
      top_SNP = SNP[is_top_snp][1L],
      top_SNP_GWS_original = any(is_top_snp & gws_original),
      top_SNP_GWS_standardized = any(is_top_snp & gws_standardized),
      n_lead_SNPs_GWS_original = sum(gws_original),
      n_lead_SNPs_GWS_standardized = sum(gws_standardized),
      any_predefined_lead_SNP_GWS_original = any(gws_original),
      any_predefined_lead_SNP_GWS_standardized = any(gws_standardized),
      best_predefined_lead_SNP_standardized = SNP[best],
      minimum_predefined_lead_SNP_P_standardized = p_standardized[best]
    )
  }, by = .(phenotype, genomic_locus, locus_id, locus_chr, locus_start, locus_end)]

  locus[, locus_retention_status := fcase(
    any_predefined_lead_SNP_GWS_original & any_predefined_lead_SNP_GWS_standardized,
    "locus_retained_at_least_one_predefined_GWS_lead_SNP",
    any_predefined_lead_SNP_GWS_original & !any_predefined_lead_SNP_GWS_standardized,
    "locus_lost_all_predefined_GWS_lead_SNPs",
    default = "locus_had_no_original_GWS_lead_SNP"
  )]

  by_phenotype <- locus[, .(
    n_source_loci = .N,
    n_top_SNPs_remained_GWS = sum(top_SNP_GWS_standardized),
    n_top_SNPs_lost_GWS = sum(top_SNP_GWS_original & !top_SNP_GWS_standardized),
    n_loci_with_any_predefined_lead_SNP_remained_GWS = sum(any_predefined_lead_SNP_GWS_standardized),
    n_loci_lost_all_predefined_GWS_lead_SNPs = sum(
      any_predefined_lead_SNP_GWS_original & !any_predefined_lead_SNP_GWS_standardized
    )
  ), by = phenotype]

  list(source = source, locus = locus, by_phenotype = by_phenotype)
}

run_effective_sample_size_standardization <- function(
    harmonized_effect_file,
    lead_manifest_file,
    sample_size_file,
    output_dir
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  effects <- fread(harmonized_effect_file, na.strings = c("NA", "", "."))
  leads <- fread(lead_manifest_file, na.strings = c("NA", "", "."))
  design <- read_sample_design(sample_size_file)
  standardized <- standardize_effects(effects, design)
  summaries <- summarize_source_loci(standardized, leads)

  phenotype_summary <- standardized[, .(
    n_unique_lead_SNPs = uniqueN(SNP),
    n_GWS_original = sum(gws_original),
    n_GWS_standardized = sum(gws_standardized)
  ), by = phenotype]
  phenotype_summary <- merge(phenotype_summary, design, by = "phenotype")

  fwrite(design, file.path(output_dir, "01_effective_sample_size_parameters.tsv"), sep = "\t")
  fwrite(standardized, file.path(output_dir, "02_standardized_results_long.tsv"), sep = "\t")
  fwrite(phenotype_summary, file.path(output_dir, "03_phenotype_summary.tsv"), sep = "\t")
  fwrite(summaries$source, file.path(output_dir, "04_source_lead_SNP_results.tsv"), sep = "\t")
  fwrite(summaries$locus, file.path(output_dir, "05_source_locus_retention.tsv"), sep = "\t")
  fwrite(summaries$by_phenotype, file.path(output_dir, "06_locus_retention_by_phenotype.tsv"), sep = "\t")
  invisible(list(standardized = standardized, locus_summary = summaries$locus))
}

run_neff_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) {
    stop("Usage: Rscript 04_standardize_effective_sample_size.R <harmonized_effects.tsv> <lead_manifest.tsv> <sample_sizes.tsv> <output_dir>")
  }
  run_effective_sample_size_standardization(args[1L], args[2L], args[3L], args[4L])
}

if (sys.nframe() == 0L) run_neff_cli()

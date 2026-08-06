#!/usr/bin/env Rscript

# Hypertension-subtype phenotype construction
#
# This public script contains only the study-specific phenotype logic. It does
# not read UK Biobank field exports, participant identifiers, withdrawal lists,
# or other individual-level source files. Run it inside the approved analysis
# environment after the general exclusions have been applied.

BP_COLUMN_MAP <- list(
  dbp = c("auto_dbp_1", "auto_dbp_2", "manual_dbp_1", "manual_dbp_2"),
  sbp = c("auto_sbp_1", "auto_sbp_2", "manual_sbp_1", "manual_sbp_2")
)

normalise_binary_flag <- function(x, column_name) {
  if (is.logical(x)) return(x)

  if (is.numeric(x) || is.integer(x)) {
    bad <- !is.na(x) & !x %in% c(0, 1)
    if (any(bad)) stop(column_name, " must contain only 0, 1, or NA.")
    return(ifelse(is.na(x), NA, x == 1))
  }

  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("1", "true", "t", "yes", "y")] <- TRUE
  out[y %in% c("0", "false", "f", "no", "n")] <- FALSE
  out[is.na(x) | y %in% c("", "na")] <- NA
  if (any(is.na(out) & !is.na(x) & nzchar(y) & y != "na")) {
    stop(column_name, " contains values that cannot be interpreted as binary.")
  }
  out
}

coerce_numeric_matrix <- function(data, columns) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop("Missing blood-pressure columns: ", paste(missing_columns, collapse = ", "))
  }
  x <- do.call(cbind, lapply(columns, function(column) data[[column]]))
  colnames(x) <- columns
  suppressWarnings(storage.mode(x) <- "double")
  x
}

row_any <- function(x) {
  apply(x, 1L, function(z) any(z, na.rm = TRUE))
}

construct_baseline_bp <- function(data, column_map = BP_COLUMN_MAP) {
  dbp <- coerce_numeric_matrix(data, column_map$dbp)
  sbp <- coerce_numeric_matrix(data, column_map$sbp)

  dbp_outlier <- row_any(dbp < 50 | dbp >= 150)
  sbp_outlier <- row_any(sbp < 70 | sbp >= 270)
  any_outlier <- dbp_outlier | sbp_outlier

  n_dbp <- rowSums(!is.na(dbp))
  n_sbp <- rowSums(!is.na(sbp))
  n_paired_readings <- rowSums(!is.na(dbp) & !is.na(sbp))

  mean_or_na <- function(x) {
    ans <- rowMeans(x, na.rm = TRUE)
    ans[is.nan(ans)] <- NA_real_
    ans
  }

  final_dbp <- mean_or_na(dbp)
  final_sbp <- mean_or_na(sbp)
  final_dbp[any_outlier] <- NA_real_
  final_sbp[any_outlier] <- NA_real_

  data.frame(
    final_sbp = final_sbp,
    final_dbp = final_dbp,
    n_sbp_readings = n_sbp,
    n_dbp_readings = n_dbp,
    n_paired_readings = n_paired_readings,
    any_bp_outlier = any_outlier,
    valid_baseline_bp = !any_outlier & !is.na(final_sbp) & !is.na(final_dbp),
    stringsAsFactors = FALSE
  )
}

define_hypertension_subtypes <- function(
    data,
    id_column = "participant_id",
    eligible_column = "preanalysis_eligible",
    medication_column = "antihypertensive_medication",
    history_column = "hypertension_history",
    age_column = "age",
    column_map = BP_COLUMN_MAP
) {
  required <- c(id_column, medication_column, history_column)
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0L) {
    stop("Missing required columns: ", paste(missing_required, collapse = ", "))
  }

  eligible <- if (eligible_column %in% names(data)) {
    normalise_binary_flag(data[[eligible_column]], eligible_column)
  } else {
    rep(TRUE, nrow(data))
  }
  medication <- normalise_binary_flag(data[[medication_column]], medication_column)
  history <- normalise_binary_flag(data[[history_column]], history_column)
  age <- if (age_column %in% names(data)) suppressWarnings(as.numeric(data[[age_column]])) else rep(NA_real_, nrow(data))

  bp <- construct_baseline_bp(data, column_map)
  untreated_and_eligible <- eligible %in% TRUE & medication %in% FALSE & bp$valid_baseline_bp

  idh <- untreated_and_eligible & bp$final_dbp >= 90 & bp$final_sbp < 140
  ish <- untreated_and_eligible & bp$final_dbp < 90 & bp$final_sbp >= 140
  sdh <- untreated_and_eligible & bp$final_dbp >= 90 & bp$final_sbp >= 140
  strict_control <- untreated_and_eligible & history %in% FALSE &
    bp$final_sbp < 120 & bp$final_dbp < 70

  subtype <- rep("not_in_case_or_reference_groups", nrow(data))
  subtype[!(eligible %in% TRUE)] <- "ineligible_before_bp_classification"
  subtype[eligible %in% TRUE & medication %in% TRUE] <- "excluded_antihypertensive_treatment"
  subtype[eligible %in% TRUE & is.na(medication)] <- "missing_treatment_status"
  subtype[eligible %in% TRUE & medication %in% FALSE & !bp$valid_baseline_bp] <- "missing_or_invalid_bp"
  subtype[strict_control] <- "CONTROL"
  subtype[idh] <- "IDH"
  subtype[ish] <- "ISH"
  subtype[sdh] <- "SDH"

  out <- data.frame(
    participant_id = data[[id_column]],
    bp,
    antihypertensive_medication = medication,
    hypertension_history = history,
    preanalysis_eligible = eligible,
    age = age,
    subtype = subtype,
    IDH_case = idh,
    ISH_case = ish,
    SDH_case = sdh,
    strict_control = strict_control,
    IDH_age_le_55_case = idh & !is.na(age) & age <= 55,
    ISH_age_gt_55_case = ish & !is.na(age) & age > 55,
    stringsAsFactors = FALSE
  )

  if (any(rowSums(out[, c("IDH_case", "ISH_case", "SDH_case", "strict_control")]) > 1L)) {
    stop("Subtype definitions are not mutually exclusive; check the input data.")
  }
  out
}

make_case_control_table <- function(classified, phenotype) {
  phenotype <- match.arg(phenotype, c("IDH", "ISH", "SDH", "IDH_age_le_55", "ISH_age_gt_55"))
  case_flag <- switch(
    phenotype,
    IDH = classified$IDH_case,
    ISH = classified$ISH_case,
    SDH = classified$SDH_case,
    IDH_age_le_55 = classified$IDH_age_le_55_case,
    ISH_age_gt_55 = classified$ISH_age_gt_55_case
  )
  control_flag <- classified$strict_control
  if (phenotype == "IDH_age_le_55") control_flag <- control_flag & classified$age <= 55
  if (phenotype == "ISH_age_gt_55") control_flag <- control_flag & classified$age > 55

  keep <- (case_flag %in% TRUE) | (control_flag %in% TRUE)
  data.frame(
    participant_id = classified$participant_id[keep],
    phenotype = as.integer(case_flag[keep]),
    stringsAsFactors = FALSE
  )
}

run_phenotype_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) {
    stop("Usage: Rscript 01_define_hypertension_subtypes.R <input.tsv> <output_dir>")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required for command-line use.")
  input <- data.table::fread(args[1L], na.strings = c("NA", "", "."), data.table = FALSE)
  classified <- define_hypertension_subtypes(input)
  dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(classified, file.path(args[2L], "baseline_subtype_classification.tsv"), sep = "\t", na = "NA")
  for (ph in c("IDH", "ISH", "SDH", "IDH_age_le_55", "ISH_age_gt_55")) {
    data.table::fwrite(
      make_case_control_table(classified, ph),
      file.path(args[2L], paste0(ph, "_case_control.tsv")),
      sep = "\t",
      na = "NA"
    )
  }
}

if (sys.nframe() == 0L) run_phenotype_cli()

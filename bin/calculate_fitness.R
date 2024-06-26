#!/usr/bin/env Rscript
#
# Purpose:
# R script to summarize sequencing read counts obtained from a CRISPRi library.
# The script summarizes count tables per sample into one main table, adds
# statistical metrics for a pairwise sample comparison using DESeq2, and calculates
# fitness scores for each gene and condition.
#
# Date: 2022-04-21
# Author: Michael Jahn, PhD
# Affilation: Science For Life Laboratory (KTH), Stockholm, Sweden

# PARSE INPUT ARGS
# ====================
#
# input arguments
args <- commandArgs(trailingOnly = TRUE)
path_samplesheet <- args[1] # file paths to sample sheet
path_counts <- args[2] # file paths to count tables
normalization <- as.logical(args[3]) # default: FALSE
gene_fitness <- as.logical(args[4]) # default: TRUE
gene_sep <- args[5] # default: "|" aka the pipe symbol
gene_controls <- args[6] # default: "" aka empty string
number_cores <- as.numeric(args[7]) # number of CPU cores


# LOAD PACKAGES
# ====================
#
library(readr)
library(tibble)
library(stringr)
library(tidyr)
library(dplyr)
library(purrr)
library(DESeq2)
library(BiocParallel)

if (normalization) {
    library(limma)
}

# DATA PREPARATION
# ====================
#
# Step 1: Load sample layout sheet - file names must be row names
df_samplesheet <- readr::read_csv(path_samplesheet, col_types = cols()) %>%
    select(all_of(c("sample", "condition", "replicate", "time", "group", "reference_group"))) %>%
    dplyr::mutate(group = factor(`group`))
stopifnot(is.numeric(df_samplesheet$time))

# Check which conditions
# - do not have differing groups/reference when zero time
# - have the minimum of two distinct time points OR
# - are compared to the zero time point of another condition OR
# - for all others skip DESeq2 and fitness calculation
df_samplesheet <- df_samplesheet %>%
    dplyr::group_by(condition) %>%
    filter(!(time == 0 & group != reference_group)) %>%
    filter(length(unique(time)) >= 2 || (time != 0 & group != reference_group)) %>%
    tibble::column_to_rownames("sample")

df_counts <- readr::read_tsv(path_counts, col_types = cols())
df_counts <- tidyr::pivot_longer(df_counts,
    cols = 3:ncol(df_counts),
    names_to = "sample", values_to = "numreads"
)

# print overview information to console
message("Number of sgRNAs detected in n samples:")
df_counts %>%
    group_by(sgRNA) %>%
    dplyr::summarize(sgRNAs_detected_in_samples = sum(numreads > 0)) %>%
    dplyr::count(sgRNAs_detected_in_samples) %>%
    dplyr::mutate(percent_total = n / sum(n) * 100) %>%
    dplyr::arrange(dplyr::desc(sgRNAs_detected_in_samples)) %>%
    print()

# input data frame must be reshaped to a 'counts matrix' with genes as rows
# and samples (conditions) as columns.
counts <- df_counts %>%
    dplyr::select(-Gene) %>%
    # spread condition over columns and sgRNAs over rows
    tidyr::pivot_wider(names_from = sample, values_from = numreads) %>%
    # remove sgRNA column, replace NA with 0
    dplyr::mutate_at(vars(-1), function(x) coalesce(x, 0)) %>%
    # add row_names from column
    tibble::column_to_rownames("sgRNA")

# DIFFERENTIAL ABUNDANCE
# ======================
#
# DESeq2 can be used to obtain fold changes and significance metrics
# for condition-wise comparisons, for details see publication:
# Love, M.I., Huber, W., Anders, S. Genome Biology, 15:550, 2014.
# (https://doi.org/10.1186/s13059-014-0550-8)

if (nrow(df_samplesheet) == 0) {
    message("No condition has sufficient groups/time points.\nFitness score calculation is omitted")
} else {
    message("Running DESeq2 for pairwise comparison.\nNote: this step can be time and computation-intense.")
    message(paste0("Number of CPU cores used for DESeq parallel execution: ", number_cores, "."))

    # 2. Meta data
    # Meta data is required to carry out the actual DESeq2 analysis
    # by 'contrasting' (comparing) selected conditions to each other.
    # We check that the order of file names corresponds to colnames of counts
    if (!all(row.names(df_samplesheet) %in% colnames(counts))) {
        stop("Not all samples listed in the sample sheet have
            corresponding read count data.")
    }
    counts <- counts[row.names(df_samplesheet)]

    # 3. Perform DESeq2 analysis
    DESeq_result <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = df_samplesheet,
        design = ~group
    ) %>%
        DESeq()

    # The syntax to call DESeq2's `results(...)` function is to use one pair of
    # contrasts `contrast("variable", "level1", "level2")`. To automate this,
    # a list of condition and reference pairs is set up from meta data
    combinations <- df_samplesheet %>%
        dplyr::select(group, reference_group) %>%
        dplyr::filter(group != reference_group) %>%
        dplyr::mutate(across(.cols = everything(), .fns = as.character)) %>%
        dplyr::distinct() %>%
        as.list() %>%
        purrr::transpose()

    # extract results for desired combinations
    DESeq_result_table <- lapply(combinations, function(l) {
        DESeq2::results(DESeq_result,
            contrast = c("group", l$group, l$reference_group),
            parallel = TRUE, BPPARAM = MulticoreParam(number_cores),
            tidy = TRUE
        ) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(group = l$group) %>%
            dplyr::rename(sgRNA = row)
    }) %>% dplyr::bind_rows()

    # MERGE DESEQ RESULTS
    # ======================
    #
    # complete metadata table with missing time zero conditions for the cases
    # where samples are used as zero time points for_multiple_different conditions
    if (!all(df_samplesheet %>% dplyr::group_by(condition) %>%
        dplyr::summarize(zero_cond = 0 %in% time) %>%
        dplyr::pull(zero_cond))
    ) {
        df_samplesheet <- dplyr::bind_rows(
            df_samplesheet,
            df_samplesheet %>% tidyr::complete(condition, time) %>%
                dplyr::filter(is.na(group), time == 0)
        )
    }

    # merge DESeq result table with meta data
    DESeq_result_table <- dplyr::select(df_samplesheet, -replicate) %>%
        dplyr::distinct() %>%
        tibble::as_tibble() %>%
        dplyr::full_join(DESeq_result_table, by = "group") %>%
        dplyr::mutate(group = as.numeric(group)) %>%
        # complete missing combinations of variables, here mostly all log2FC
        # values (0) for the reference conditions
        tidyr::complete(sgRNA, tidyr::nesting(condition, time, group, reference_group)) %>%
        dplyr::mutate(
            log2FoldChange = replace_na(log2FoldChange, 0),
            lfcSE = replace_na(lfcSE, 0),
            pvalue = replace_na(pvalue, 1),
            padj = replace_na(padj, 1)
        ) %>%
        dplyr::filter(!is.na(sgRNA))

    # NORMALIZATION
    # ======================
    #
    # input data frame must be reshaped to a 'counts matrix' with genes as rows
    # and samples (conditions) as columns.
    # in order to do this, construct a normalization function that takes three
    # columns as input, the numeric variable to be normalized, the conditioning variable
    # (character or factor), and an ID that identifies each observation
    # (barcode/mutant/sgRNA)
    if (normalization) {
        message("Performing quantile normalization on per sample read counts.")

        apply_norm <- function(id, cond, var) {
            df_orig <- tibble(id = id, cond = cond, var = var)
            df_new <- pivot_wider(df_orig, names_from = cond, values_from = var) %>%
                column_to_rownames("id") %>%
                as.matrix() %>%
                limma::normalizeBetweenArrays(method = "quantile") %>%
                as_tibble(rownames = "id") %>%
                pivot_longer(-id, names_to = "cond", values_to = "var_norm")
            left_join(df_orig, df_new, by = c("id", "cond")) %>% pull(var_norm)
        }

        # apply normalization
        DESeq_result_table <- DESeq_result_table %>%
            mutate(FoldChange = 2^log2FoldChange) %>%
            group_by(time) %>%
            mutate(
                FoldChange_norm = apply_norm(sgRNA, condition, FoldChange),
                log2FoldChange = log2(FoldChange_norm)
            ) %>%
            ungroup() %>%
            select(-FoldChange, -FoldChange_norm)
    }

    # CALCULATE FITNESS SCORE
    # =======================
    #
    # Here we define fitness score as the area under/over the curve for log2 fold change
    # over time. Enrichment will result in a positive score, depletion
    # in a negative score. The fitness score is normalized to the maximum time
    # for a particular condition, and is therefore independent of the duration
    # of the cultivations. Requires at least 2 time points
    calc_auc <- function(x, y) {
        sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
    }
    message("Calculating sgRNA fitness score.")
    DESeq_result_table <- DESeq_result_table %>%
        dplyr::arrange(sgRNA, condition, time) %>%
        dplyr::group_by(sgRNA, condition) %>%
        dplyr::mutate(fitness = calc_auc(time, log2FoldChange) / (max(time) / 2))


    # CALCULATE SGRNA CORRELATION AND EFFICIENCY
    # ==========================================
    #
    # Different sgRNAs per gene can have different repression efficiency.
    # To assess sgRNA quality, two metrics are added to the main table,
    # A) sgRNA correlation = Pearson correlation coeff. of each sgRNA with the others of same gene.
    #        A score between 0 and 1.
    # B) sgRNA efficiency = median absolute fitness of an sgRNA over all observations [conditions],
    #        divided by maximum fitness of an sgRNA. A score between 0 and 1.
    if (as.logical(gene_fitness)) {
        message("Calculating sgRNA efficiency and correlation.")
        determine_corr <- function(index, value, condition, time) {
            # make correlation matrix
            df <- data.frame(index = index, value = value, condition = condition, time = time)
            cor_matrix <- tidyr::pivot_wider(df, names_from = c("condition", "time"), values_from = value) %>%
                dplyr::arrange(index) %>%
                tibble::column_to_rownames("index") %>%
                as.matrix() %>%
                t() %>%
                stats::cor(method = "pearson")
            # determine weights
            weights <- cor_matrix %>%
                replace(., . == 1, NA) %>%
                apply(2, function(x) median(x, na.rm = TRUE)) %>%
                sapply(function(x) {
                    (x + 1) / 2
                }) %>%
                tibble::enframe("index", "weight") %>%
                dplyr::mutate(index = as.numeric(index)) %>%
                dplyr::mutate(weight = replace(weight, is.na(weight), 1))
            # return vector of weights the same order and length
            # as sgRNA index vector
            dplyr::left_join(df, weights, by = "index") %>% pull(weight)
        }

        DESeq_result_table <- DESeq_result_table %>%
            # split sgRNA names into target gene and position
            tidyr::separate(sgRNA,
                into = c("sgRNA_target", "sgRNA_position"), sep = paste0("\\", gene_sep),
                remove = FALSE
            ) %>%
            dplyr::group_by(sgRNA_target) %>%
            dplyr::mutate(
                sgRNA_position = as.numeric(sgRNA_position),
                sgRNA_index = sgRNA_position %>% as.factor() %>% as.numeric(),
                sgRNA_correlation = determine_corr(sgRNA_index, log2FoldChange, condition, time),
                sgRNA_efficiency = ave(fitness, sgRNA_position, FUN = function(x) median(abs(x))) %>%
                    {
                        . / max(.)
                    }
            )

        # SUMMARIZE SGRNA FITNESS TO GENE FITNESS
        # =======================================
        #
        # calculate gene fitness as weighted mean of sgRNA fitness, see README.md
        # for details and exact formula
        message("Calculating gene fitness and gene log2 fold change.")
        if (gene_controls != "") {
            df_controls <- DESeq_result_table %>%
                ungroup() %>%
                filter(str_detect(sgRNA_target, gene_controls))
        }
        get_controls <- function(condition) {
            if (gene_controls != "") {
                filter(df_controls, condition == condition)$fitness
            } else {
                NULL
            }
        }

        DESeq_result_table <- dplyr::left_join(
            DESeq_result_table,
            DESeq_result_table %>%
                dplyr::group_by(sgRNA_target, condition, time) %>%
                dplyr::summarize(
                    .groups = "keep",
                    # gene log2 FC
                    wmean_log2FoldChange = weighted.mean(log2FoldChange, sgRNA_correlation * sgRNA_efficiency),
                    sd_log2FoldChange = sd(log2FoldChange),
                    # gene fitness
                    wmean_fitness = weighted.mean(fitness, sgRNA_correlation * sgRNA_efficiency),
                    sd_fitness = sd(fitness),
                    # apply Wilcoxon rank sum test against Null hypothesis (fitness ~ 0) or
                    # user-supplied set of controls
                    p_fitness = stats::wilcox.test(
                        x = fitness,
                        y = get_controls(condition),
                        alternative = "two.sided"
                    )$p.value
                ),
            by = c("sgRNA_target", "condition", "time")
        ) %>%
            group_by(condition, time) %>%
            mutate(
                p_fitness_adj = stats::p.adjust(p_fitness, method = "BH"),
                comb_score = abs(wmean_fitness) * -log10(p_fitness_adj)
            )
    }
}

# EXPORT PROCESSED DATA
# =====================
#
# Save result tables to output folder, in Rdata format
if (nrow(df_samplesheet) > 0) {
    message("Saving 'result.Rdata'")
    save(DESeq_result_table, file = "result.Rdata")
    message("Saving 'result.tsv'")
    readr::write_tsv(DESeq_result_table, file = "result.tsv")
}

#!/usr/bin/env Rscript

# Purpose:
# R script to import counts files from subread/featureCounts module;
# All counts files are then merged to one master count table as
# input for Mageck. This script is part of the nextflow pipeline
# 'nf-core-crispriscreen' for processing of CRISPRi library data.
# The pipeline can be found at https://github.com/MPUSP/nf-core-crispriscreen.
#
# Date: 2023-02-14
# Author: Michael Jahn, PhD
# Affilation: Science For Life Laboratory (KTH), Stockholm, Sweden

# input parameters
args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1] # file path to sample sheet
path_counts <- args[2] # file paths to count tables
gene_sep <- args[3] # default: "|" aka the pipe symbol
input_design <- args[4] # file path to custom design matrix, default: null


### PREPARE MASTER COUNT TABLE

# check gene guide separator
if (gene_sep %in% c("|", "_", ":", ";", ".", "/", " ")) {
    gene_sep <- paste0("\\", gene_sep)
}

# read count tables to list
list_files <- list.files(
    path = getwd(),
    full.names = TRUE,
    pattern = ".featureCounts.txt$"
)

list_counts <- lapply(list_files, function(x) {
    df <- read.delim(x, skip = 1, stringsAsFactors = FALSE)
    df[c(1, 7)]
})

# combine all count tables and remove duplicated annotation
df_combined <- do.call(cbind, list_counts)
stopifnot(df_combined[[1]] == df_combined[[3]])
df_combined <- df_combined[!duplicated(colnames(df_combined))]

# separate sgRNA and gene names, reorder columns
sgrna <- df_combined[[1]]
df_combined$Gene <- substr(sgrna, 1, as.numeric(regexpr(gene_sep, sgrna) - 1))
colnames(df_combined)[1] <- "sgRNA"
df_combined <- df_combined[c(1, ncol(df_combined), 2:(ncol(df_combined) - 1))]
colnames(df_combined)[3:ncol(df_combined)] <- gsub(
    "(\\_T[0-9]+)?\\.bam$", "",
    colnames(df_combined)[3:ncol(df_combined)]
)

write.table(
    x = df_combined, file = "all_counts.tsv", sep = "\t",
    row.names = FALSE, col.names = TRUE, quote = FALSE
)

### PREPARE DESIGN MATRIX

if (input_design == "\\|") {
    df_sample <- read.csv(samplesheet, stringsAsFactors = FALSE)
    df_design <- data.frame(
        sample = df_sample[["sample"]],
        baseline = 1,
        time = df_sample[["time"]],
        common = as.numeric(df_sample[["time"]] != 0)
    )

    # add 1 column for each condition/treatment with
    # binary encoding of samples that determine condition
    for (cond in unique(df_sample$condition)) {
        df_design[[cond]] <- as.numeric(
            (df_sample$condition %in% cond) & (df_sample$time != 0)
        )
    }
} else {
    df_design <- read.delim(input_design,
        stringsAsFactors = FALSE,
        header = TRUE,
        row.names = NULL
    )
}

# test that design matrix corresponds to counts table
stopifnot(all(df_design$sample %in% colnames(df_combined)[-c(1, 2)]))

# export design matrix
write.table(
    x = df_design, file = "design.tsv", sep = "\t",
    row.names = FALSE, col.names = TRUE, quote = FALSE
)

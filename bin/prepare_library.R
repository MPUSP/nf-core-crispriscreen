#!/usr/bin/env Rscript

# Purpose:
# R script to import a *.fasta file and convert it to an *.saf genome
# annotation file. This script is part of the nextflow pipeline
# 'nf-core-crispriscreen' for processing of CRISPRi library data.
# The pipeline can be found at https://github.com/m-jahn/nf-core-crispriscreen.
# It uses only base R to keep dependencies low.
#
# Date: 2022-04-11
# Author: Michael Jahn, PhD
# Affilation: Science For Life Laboratory (KTH), Stockholm, Sweden

# input parameter: path to fasta file
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]

# read and process fasta file
fasta_df <- read.delim(fasta_file, header = FALSE)
saf_df <- data.frame(
    GeneID = fasta_df[seq(1, nrow(fasta_df), 2), 1],
    Sequence = fasta_df[seq(2, nrow(fasta_df), 2), 1]
)

# check for duplications
stopifnot(!any(duplicated(fasta_df$GeneID)))
stopifnot(!any(duplicated(fasta_df$Sequence)))

# add remaining columns
saf_df$GeneID <- gsub(">", "", saf_df$GeneID)
saf_df$Chr <- saf_df$GeneID
saf_df$Start <- 1
saf_df$End <- sapply(saf_df$Sequence, nchar)
saf_df$Strand <- "*"
saf_df <- saf_df[c("GeneID", "Chr", "Start", "End", "Strand", "Sequence")]

saf_name <- gsub("fasta$", "saf", basename(fasta_file))
write.table(x = saf_df, file = saf_name, sep = "\t", row.names = FALSE,
    col.names = FALSE, quote = FALSE)

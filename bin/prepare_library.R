#!/usr/bin/env Rscript

# Purpose:
# R script to import a *.fasta file and convert it to an *.saf genome
# annotation file. This script is part of the nextflow pipeline
# 'nf-core-crispriscreen' for processing of CRISPRi library data.
# The pipeline can be found at https://github.com/MPUSP/nf-core-crispriscreen.
#
# Date: 2022-04-11
# Author: Michael Jahn, PhD
# Affilation: Science For Life Laboratory (KTH), Stockholm, Sweden

# input parameters
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1] # path to fasta file, mandatory
gene_controls <- args[2] # pattern for control barcodes, default: "" aka empty string

# read and process fasta file
fasta_df <- read.delim(fasta_file, header = FALSE)
## identify lines which start with ">" - if there are none, raise error (then probably no .fasta file!) - catch cases in which nucleotide sequence is spread over several lines
number_headers <- sum(grepl(">", fasta_df[,1]))
stopifnot(!number_headers==0)
header_vec <- 1:number_headers
nt_seq_vec <- 1:number_headers
tmp_nt_seq <- ""
counter_headers <- 1
for(line_i in 1:nrow(fasta_df)){
  if(grepl(">", fasta_df[line_i,])){
    if(counter_headers == 1){
      header_vec[counter_headers] <- fasta_df[line_i,]
      counter_headers = counter_headers + 1
    } else {
      nt_seq_vec[counter_headers-1] <- tmp_nt_seq
      header_vec[counter_headers] <- fasta_df[line_i,]
      counter_headers = counter_headers + 1
    }
    tmp_nt_seq = ""
  } else {
    tmp_nt_seq <- paste(tmp_nt_seq, trimws(fasta_df[line_i,]), sep="")
  }
  if(line_i == nrow(fasta_df)){
    nt_seq_vec[counter_headers-1] <- tmp_nt_seq
  }
}

saf_df <- data.frame(
    GeneID = header_vec,
    Sequence = nt_seq_vec
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

# replace file ending by .saf (everything after last . in file name is interpreted as file ending) - if there is no file ending, append ".saf"
if (length(strsplit(fasta_file, "\\.")[[1]][-1]) > 0) {
  saf_name <- gsub(strsplit(fasta_file, "\\.")[[1]][-1], "saf", basename(fasta_file))
  } else {
  saf_name <- paste(fasta_file, ".saf", sep="")
}
write.table(
    x = saf_df, file = saf_name, sep = "\t", row.names = FALSE,
    col.names = FALSE, quote = FALSE
)

# optionally identify controls and save them for reference
if (gene_controls != "") {
    ctrl_hits <- grepl(gene_controls, saf_df$GeneID)
    if (!sum(ctrl_hits)) {
        stop(paste0("No barcode name matches supplied pattern '", gene_controls, "'."))
    }
    saf_df_ctrl <- subset(saf_df, ctrl_hits)
    ctrl_name <- gsub(".fasta$", "_controls.tsv", basename(fasta_file))
    write.table(
        x = saf_df_ctrl, file = ctrl_name, sep = "\t", row.names = FALSE,
        col.names = FALSE, quote = FALSE
    )
}


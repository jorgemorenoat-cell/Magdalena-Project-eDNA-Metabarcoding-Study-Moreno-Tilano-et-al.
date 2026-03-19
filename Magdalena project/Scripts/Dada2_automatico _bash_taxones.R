#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(dplyr)
})

# ------------------------------
# Argumentos desde bash
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]   # Carpeta con fastq.gz
marker <- args[2]      # Ej: "teleo" o "vert01"

## Bases de referencia (ajusta rutas según corresponda)
#teleo
#ref_tax_all <- ("../../../../../dada2_pipeline/utils/Dada_mitofish_teleo_taxonomy")
#vert01
ref_tax_all <- ("../../../../dada2_pipeline/utils/Dada_mitofish_midori_mitofish_vert_taxonomy")

# Nombre del subset
subset_name <- basename(normalizePath(input_dir))

# Carpeta de salida
output_dir <- file.path("outputs", "Taxa_final", subset_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------
# Leer fastq y pipeline DADA2
# ------------------------------
fnFs <- sort(list.files(input_dir, pattern=".R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern=".R2.fastq.gz", full.names = TRUE))

filtpath <- file.path(tempdir(), paste0("filt_", subset_name))
dir.create(filtpath, showWarnings = FALSE)
filtFs <- file.path(filtpath, basename(fnFs))
filtRs <- file.path(filtpath, basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxEE=c(2,2), truncQ=11, maxN=0, rm.phix=TRUE,
                     compress=FALSE, verbose=TRUE, multithread=TRUE)

set.seed(100)
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "\\."), `[`, 1)
names(filtFs) <- sample.names
names(filtRs) <- sample.names

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  mergers[[sam]] <- mergePairs(ddF, derepF, ddR, derepR, minOverlap=20)
}

seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

# ------------------------------
# Asignación taxonómica
# ------------------------------
taxa <- assignTaxonomy(seqtab_nochim, ref_tax_all, multithread = TRUE, tryRC = TRUE)
taxa_df <- as.data.frame(taxa)

# ------------------------------
# Agrupar por taxonomía y contar
# ------------------------------
taxa_consolidados <- taxa_df %>%
  distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all = FALSE)

n_taxa <- nrow(taxa_consolidados)

# ------------------------------
# Guardar resultados
# ------------------------------
# CSV de taxones (concatenando ASVs que caen en el mismo taxón)
csv_out <- file.path(output_dir, paste0("Taxa_", subset_name, ".csv"))
write.csv(taxa_consolidados, file = csv_out, row.names = FALSE)

# CSV acumulado con número de taxones
counts_out <- file.path("outputs", "Taxa_final", "Taxa_counts.csv")
if (!file.exists(counts_out)) {
  cat("Subset,Num_Taxa\n", file=counts_out)
}
cat(paste0(subset_name, ",", n_taxa, "\n"), file=counts_out, append=TRUE)

cat("Pipeline completado para", subset_name, ":", n_taxa, "taxones\n")

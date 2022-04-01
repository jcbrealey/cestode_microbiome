#!/usr/bin/Rscript

### Rscript for running DADA2 pipeline with Rv4.0.x
### Big data approach
### Step 2: infer sequence variants
library(dada2)
library(phyloseq)
print("Packages: dada2, phyloseq; Versions:")
packageVersion("dada2"); packageVersion("phyloseq")

# Set up directories and file names
filt.dir <- "/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/Pd2_dadaproc_200814"
out.dir <- "/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/Pd3_dadaASV_200814"
filtFs <- list.files(filt.dir, pattern="F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filt.dir, pattern="R_filt.fastq.gz", full.names = TRUE)
tax_file <- "/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/SILVA_138_SSU_dada2/silva_nr_v138_train_set_Limborg.fa.gz"

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

dir.create(file.path(out.dir))

set.seed(923)
# Learn error rates - use subset of nt to speed things up (see tutorial)
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
# 14Aug2020 v2: try maxMismatch=2 instead of 0
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap=12, maxMismatch=2)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nc <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

print("raw ASV table dimensions")
dim(seqtab)
print("ASV table post-chimera removal")
dim(seqtab.nc)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nc, tax_file, multithread=TRUE)

# sanity check
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) #default ordered by abundance in entire dataset

# Save DADA2 outputs
saveRDS(seqtab, paste0(out.dir,"/TW16S_dada2_seqtab_200814_v2.rds"))
saveRDS(seqtab.nc, paste0(out.dir,"/TW16S_dada2_seqtabnc_200814_v2.rds"))
saveRDS(taxa, paste0(out.dir,"/TW16S_dada2_seqtabnc_TAXA_200814_v2.rds"))

# get unknown counts
ps <- phyloseq(otu_table(seqtab.nc, taxa_are_rows=FALSE), tax_table(taxa))
ps.kingdom <- tax_glom(ps, taxrank="Kingdom", NArm=FALSE)
df.king <- as.data.frame(otu_table(ps.kingdom))
names(df.king) <- as.data.frame(tax_table(ps.kingdom))$Kingdom
df.king$total.ct <- rowSums(df.king)

# original counts
raw.cts <- read.csv(paste0(filt.dir,"/log_dada2_filt_200814_maxEE3.csv"), row.names=1)
row.names(raw.cts) <- sapply(strsplit(row.names(raw.cts), "_"), `[`, 1) 
raw.cts$SampleID <- row.names(raw.cts)
names(raw.cts) <- c("Preprocessed","Filtered","SampleID")

df.king$SampleID <- row.names(df.king)
df.king <- merge(df.king, raw.cts, by = "SampleID")
## code from here failed, not sure why
df.king$Assigned.per <- (df.king$NA. - df.king$total.ct) / df.king$Preprocessed * 100
df.king$Run <- "Aug14_v2"

# save count ouptput
write.csv(df.king, paste0(out.dir,"log_dada2_inferASV_200814_maxEE3_mm2.csv"), row.names=T, quote=F)

save.image(paste0(out.dir,"Rworkspace_dada2_inferASV_200814_maxEE_mm2.RData"))

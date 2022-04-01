#!/usr/bin/Rscript

### Rscript for running DADA2 pipeline with Rv4.0.x
### Big data approach
### Step 3: assign taxonomy using custom/edited database
library(dada2)
library(phyloseq)
print("Packages: dada2, phyloseq; Versions:")
packageVersion("dada2"); packageVersion("phyloseq")

# Using previously inferred ASV table (200814_v2)
# but assigned taxonomy with an updated database including 
# tapeworm and salmon SSU sequences downloaded from SILVA v138

# Set up directories and file names
out.dir <- "/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/Pd3_dadaASV_200814"
seqtab_file <- paste0(out.dir,"/TW16S_dada2_seqtabnc_200814_v2.rds")
tax_file1 <- "reference_databases/silva_nr_v138_custom_200821_train_set.fa.gz"
tax_file2 <- "reference_databases/silva_nr_v138_custom_200821_species_assignment.fa.gz"

set.seed(923)

# Read in seqtab file
seqtab.nc <- readRDS(seqtab_file)
print("ASV table post-chimera removal")
dim(seqtab.nc)

# Assign taxonomy to genus-level
taxa <- assignTaxonomy(seqtab.nc, tax_file1, multithread=TRUE)

# Assign species-level taxonomy
taxa2 <- addSpecies(taxa, tax_file2)

# sanity check
taxa.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) #default ordered by abundance in entire dataset

# Save taxa outputs
saveRDS(taxa, paste0(out.dir,"/TW16S_dada2_seqtabnc_TAXA_200821_genus.rds"))
saveRDS(taxa2, paste0(out.dir,"/TW16S_dada2_seqtabnc_TAXA_200821_species.rds"))


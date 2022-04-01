#!/usr/bin/Rscript

### Rscript for running DADA2 pipeline with Rv4.0.x
### Big data approach
### Step 1: filtering
library(dada2)
packageVersion("dada2")

# Set up directories and file names
path.top <- "/groups/hologenomics/jaelleb/TAPEWORM_16S_200729"
in.dir <- paste0(path.top,"/P1_adrm_200811_PE")
out.dir <- paste0(path.top,"/Pd2_dadaproc_200814")
fastqFs <- sort(list.files(in.dir, pattern="adrm.pair1.truncated"))
fastqRs <- sort(list.files(in.dir, pattern="adrm.pair2.truncated"))

sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filt.out <- filterAndTrim(fwd=file.path(in.dir, fastqFs),
              filt=file.path(out.dir, paste0(sample.names,"_F_filt.fastq.gz")),
              rev=file.path(in.dir, fastqRs), 
              filt.rev=file.path(out.dir, paste0(sample.names,"_R_filt.fastq.gz")),
              truncLen=0, maxEE=c(3,3), truncQ=2, maxN=0, trimLeft=c(17,20), minLen=200,
              rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)

write.csv(filt.out,paste0(out.dir,"/log_dada2_filt_200814_maxEE3.csv"), row.names=T, quote=F)
.libPaths("C:/Users/jaelleb/R/Rv4.1.2_libs")
library(ggplot2)
library(reshape2)
library(phyloseq)
library(decontam)
library(vegan)
library(RColorBrewer)
library(microbiome)
library(cowplot)
library(Maaslin2)
library(MicEco)
library(hillR)
library(ggVennDiagram)
set.seed(835)

####### HoloFish Tapeworm 16S DADA2 + phyloseq analyses #######

setwd("C:/Users/jaelleb/WORK/R_analyses_Cdrive/RCdr_Tapeworm/cestode_microbiome")

#### Useful functions ####
rank_summary <- function(ps, rank, NArm=FALSE){
  # summarises at specified rank
  # takes phyloseq object, returns dataframe
  ps.rank <- tax_glom(ps, taxrank = rank, NArm = NArm)
  df.rank <- as.data.frame(otu_table(ps.rank))
  names(df.rank) <- as.data.frame(tax_table(ps.rank))[,rank]
  names(df.rank)[is.na(names(df.rank))] <- "Unknown"
  df.rank$Total.ct <- rowSums(df.rank)
  return(df.rank)
}

compare_pcr_reps <- function(seqtab, samples, ntimes = 1, nreps = 2){
  # compares ASV detection in PCR replicates
  # keeps ASVs present in nreps more than ntimes
  for (s in samples){
    rowids <- grep(s,rownames(seqtab))
    seqtab[rowids, names(which(colSums(seqtab[rowids,] > ntimes) < nreps))] <- 0
  }
  seqtab <- seqtab[,colSums(seqtab) > 0]
  return(seqtab)
}

taxonomy_rank <- function(tax.row){
  # returns the taxonomic rank at which an ASV was assigned
  ranks <- data.frame(Level=seq(1,7,1),
                      Rank=c("Unassigned","Kingdom","Phylum","Class","Order","Family","Genus"))
  
  l <- which(is.na(tax.row))
  if(length(l) == 0){
    r <- "Species"
  } else {
    r <- ranks[which(ranks$Level == l[1]),"Rank"]
  }
  return(as.character(r))
}

top_taxa <- function(ct.df, sample_var = "SampleID.R", top = 9) {
  # calculates the top x taxa for ASVs summarised at a taxonomic level
  # ct.df = dataframe: samples ~ sum_rank, value = count
  sums <- colSums(ct.df[,which(!names(ct.df) %in% c(sample_var,"Unassigned"))])
  sums.desc <- names(sums[order(sums, decreasing = T)])
  sums.1 <- ct.df[,c(sample_var,sums.desc[1:top])]
  sums.2 <- data.frame(Other=rowSums(ct.df[,sums.desc[(top+1):length(sums.desc)]]))
  if("Unassigned" %in% names(ct.df)){
    top.df <- cbind.data.frame(sums.1, sums.2, Unassigned=ct.df$Unassigned)
  } else {
    top.df <- cbind.data.frame(sums.1, sums.2)
  }
  return(top.df)
}

tab_2_df <- function(seqtab, meta.df, relab = TRUE, merge.meta = TRUE, sID){
  # take seqtab matrix and return dataframe, optionally with metadata
  # seqtab = matrix, rownames = SampleID, colnames = ASV.ID
  if(relab){
    df <- prop.table(seqtab, margin = 1)
  }
  df <- as.data.frame(df)
  df$SampleID <- rownames(seqtab)
  names(df)[ncol(df)] <- sID
  if(merge.meta){
    df <- merge(df, meta, by = sID)
  }
  return(df)
}
print_fasta <- function(seqids, asv.seqs){
  for (i in seqids) {
    n <- paste(">",i, sep = "")
    s <- as.character(asv.seqs[which(asv.seqs$ASV.ID == i),"Sequence"])
    writeLines(c(n,s))
  }
}
# emulate ggplot2 default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#### Set up data ####
# Supp table S1
meta.holofish <- read.csv("data/HoloFish_metadata_individuals.csv")
# Supp table S2
meta <- read.csv("data/HoloFish_metadata_samples_cestode.csv")

row.names(meta) <- meta$SampleID.R
meta$Cestode.index <- factor(meta$Cestode.index, levels = c(0,1,2,3), ordered = T)
meta$Sample.type <- factor(meta$Sample.type, levels = c("salmon_gut","cestode_wash","cestode_body","c_mock","c_blank","c_ddH2O"))

row.names(meta.holofish) <- meta.holofish$Fish.ID
meta.holofish$Cestode.index <- factor(meta.holofish$Cestode.index, levels = c(0,1,2,3), ordered = T)

samples <- unique(meta$SampleID)
samples.reps <- unique(meta$SampleID.R)

mock.theory <- read.csv("data/mock_community_composition.csv", sep = ",", skip = 1, header = T)
# exclude fungi from mock that won't be picked up by 16S
mock.theory <- mock.theory[which(!mock.theory$Genus %in% c("Saccharomyces","Cryptococcus")),]

asv.tab <- readRDS("data/ASV_table_nochimeras.rds")
asv.seqs <- data.frame(ASV.ID=paste("seq",seq(1,ncol(asv.tab),1), sep=""), 
                       Sequence=colnames(asv.tab))
colnames(asv.tab) <- as.character(asv.seqs$ASV.ID)

asv.taxa <- readRDS("data/ASV_taxa_nochimeras_species.rds")
rownames(asv.taxa) <- sapply(rownames(asv.taxa), function(x){
  asv.seqs[which(asv.seqs$Sequence == x), "ASV.ID"]
})
asv.taxa.df <- as.data.frame(asv.taxa)
asv.taxa.df$ASV.ID <- row.names(asv.taxa.df)
asv.taxa.df$rank <- sapply(seq(1,nrow(asv.taxa.df),1), function(i){
  taxonomy_rank(asv.taxa.df[i,])
})
asv.taxa.df$sum_Genus <- sapply(seq(1,nrow(asv.taxa.df),1), function(i){
  if(asv.taxa.df[i,"rank"] %in% c("Genus","Species")){
    as.character(asv.taxa.df[i,"Genus"])
  } else if(asv.taxa.df[i,"rank"] == "Unassigned") {
    "Unassigned"
  } else {
    l <- substring(as.character(asv.taxa.df[i,"rank"]),1,1)
    t <- as.character(asv.taxa.df[i,as.character(asv.taxa.df[i,"rank"])])
    paste(l,"__",t, sep = "")
  }
})

ps <- phyloseq(otu_table(asv.tab, taxa_are_rows=FALSE),
               tax_table(asv.taxa),
               sample_data(meta)
)

#### QC: Mock community ####
# Compare actual composition with theoretical composition
unfilt.genus.sum <- as.data.frame(t(asv.tab))
unfilt.genus.sum$ASV.ID <- row.names(unfilt.genus.sum)
unfilt.genus.sum <- merge(unfilt.genus.sum, asv.taxa.df[,c("ASV.ID","sum_Genus")], by = "ASV.ID")
unfilt.genus.sum <- melt(unfilt.genus.sum, id.vars = "sum_Genus", measure.vars = rownames(asv.tab), 
                         variable.name = "SampleID.R", value.name = "count")
unfilt.genus.sum <- dcast(unfilt.genus.sum, SampleID.R ~ sum_Genus, value.var = "count", fun.aggregate = sum)

mock.comp <- top_taxa(unfilt.genus.sum[grep("Mock",unfilt.genus.sum$SampleID.R),], top = 8)
names(mock.comp)[3] <- "Escherichia"
mock.comp <- rbind.data.frame(mock.comp[,which(names(mock.comp) != "Unassigned")], data.frame(
  SampleID.R="Mock_theory",
  Enterococcus=mock.theory[mock.theory$Genus == "Enterococcus", "rRNA_16S"]/100,
  Escherichia=mock.theory[mock.theory$Genus == "Escherichia", "rRNA_16S"]/100,
  Salmonella=mock.theory[mock.theory$Genus == "Salmonella", "rRNA_16S"]/100,
  Staphylococcus=mock.theory[mock.theory$Genus == "Staphylococcus", "rRNA_16S"]/100,
  Bacillus=mock.theory[mock.theory$Genus == "Bacillus", "rRNA_16S"]/100,
  Pseudomonas=mock.theory[mock.theory$Genus == "Pseudomonas", "rRNA_16S"]/100,
  Lactobacillus=mock.theory[mock.theory$Genus == "Lactobacillus", "rRNA_16S"]/100,
  Listeria=mock.theory[mock.theory$Genus == "Listeria", "rRNA_16S"]/100,
  Other=0.00)
)

# Supp figure: mock composition
ggplot(melt(mock.comp, id.vars = "SampleID.R", measure.vars = c(2:ncol(mock.comp)),
            variable.name = "Genus", value.name = "rRNA_16S"), 
       aes(SampleID.R, rRNA_16S, fill = Genus))+
  geom_col(position="fill", width = 0.9)+
  scale_fill_manual(values=c(brewer.pal(8, "Spectral"), "#969696"))+
  scale_x_discrete(
    labels = c("Mock1 \n1x","Mock1R \n1x","Mock2 \n1x","Mock2R \n1x",
               "Mock3 \n2x","Mock3R \n2x","Mock4 \n3x","Mock4R \n3x",
               "Mock_theory")
  )+
  geom_vline(xintercept = 8.5, size = 1)+
  labs(x = "Sample", y = "Relative abundance")+
  theme_classic()+
  theme(axis.text = element_text(size = 11, colour = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size = 11, colour = "black"), legend.title = element_text(size = 12))
if(!file.exists("figures/Supp_additional_fig_mock_composition.png")){
  ggsave("figures/Supp_additional_fig_mock_composition.png", units = "in", dpi = 600, width = 9.5, height = 4.5)
}

#### QC: PCR false positives identification ####
asv.tab.filt <- compare_pcr_reps(asv.tab, samples, ntimes = 1, nreps = 2)
pcr.false.pos <- setdiff(colnames(asv.tab),colnames(asv.tab.filt))

#### QC: contaminant identification ####
contam.freq <- isContaminant(ps, method="frequency", conc="Library.PCR.concentration")
table(contam.freq$contaminant)
#FALSE  TRUE 
# 1160    12
asv.taxa.df[row.names(contam.freq[which(contam.freq$contaminant),]),
            c("Genus","Species")]
contams.freq <- row.names(contam.freq[which(contam.freq$contaminant),])

contam.prev <- isContaminant(ps, method="prevalence", neg="Lab.neg")
table(contam.prev$contaminant)
#"FALSE  TRUE 
# 1152    20 "
asv.taxa.df[row.names(contam.prev[which(contam.prev$contaminant),]),
            c("Genus","Species")]
# many ASVs from mock

ps.unfilt.mock <- subset_samples(ps, Sample.type=="c_mock")
ps.unfilt.mock <- filter_taxa(ps.unfilt.mock, function(x) sum(x) > 0, TRUE)
ps.unfilt.mock.rel <- transform_sample_counts(ps.unfilt.mock, function(OTU) OTU/sum(OTU))
# 29 taxa and 8 samples

taxa.mockseqs <- row.names(tax_table(ps.unfilt.mock))
taxa.mockseqs <- taxa.mockseqs[which(taxa.mockseqs %in% asv.taxa.df[which(asv.taxa.df$Genus %in% c(mock.theory$Genus,"Escherichia/Shigella")),"ASV.ID"])]

# Use mock samples to identify probable cross-contamination into other samples
contam.prev.ids <- row.names(contam.prev[which(contam.prev$contaminant),])
contams.blanks <- contam.prev.ids[which(!contam.prev.ids %in% taxa.mockseqs)]

#### QC: Filter contaminants etc out of all samples ####
asv.tab.filt <- asv.tab.filt[,which(!colnames(asv.tab.filt) %in% unique(c(contams.freq, contams.blanks)))]

## set mock sequences to zero in all samples
asv.tab.filt[grep("Mock",rownames(asv.tab.filt), invert = T),
             which(colnames(asv.tab.filt) %in% taxa.mockseqs)] <- 0

asv.tab.filt <- asv.tab.filt[which(rowSums(asv.tab.filt) > 0),]
# both replications for blank1 and blank3 now empty

## sanity check: compare abundance of contams etc
contam.sums <- data.frame(
  SampleID.R=rownames(asv.tab),
  PCR.FP=rowSums(prop.table(asv.tab, margin = 1)[,pcr.false.pos]),
  contams=rowSums(prop.table(asv.tab, margin = 1)[,unique(c(contams.freq, contams.blanks))]),
  mockcross=rowSums(prop.table(asv.tab, margin = 1)[,taxa.mockseqs])
)
contam.sums <- merge(contam.sums, meta, by = "SampleID.R")
ggplot(contam.sums, aes(Sample.type, PCR.FP))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Sample.type), size = 3, position = position_jitterdodge())+
  scale_y_continuous(limits = c(0,1))+
  scale_x_discrete(labels = c("salmon gut","cestode wash","cestode body","mock","blanks","water"))+
  labs(title = "PCR false positives", x = "Sample type", y = "Relative abundance")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        legend.position = "none")
ggplot(contam.sums, aes(Sample.type, contams))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Sample.type), size = 3, position = position_jitterdodge())+
  scale_y_continuous(limits = c(0,1))+
  scale_x_discrete(labels = c("salmon gut","cestode wash","cestode body","mock","blanks","water"))+
  labs(title = "True contaminants", x = "Sample type", y = "Relative abundance")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        legend.position = "none")
ggplot(contam.sums, aes(Sample.type, mockcross))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Sample.type), size = 3, position = position_jitterdodge())+
  scale_y_continuous(limits = c(0,1))+
  scale_x_discrete(labels = c("salmon gut","cestode wash","cestode body","mock","blanks","water"))+
  labs(title = "Cross contamination from mock", x = "Sample type", y = "Relative abundance")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
        legend.position = "none")

## filter out Eukaryote, Mitochondria and Chloroplast sequences
excluded.seqs <- unique(c(
  asv.taxa.df[which(asv.taxa.df$ASV.ID %in% colnames(asv.tab.filt) & 
                      asv.taxa.df$Family == "Mitochondria"),"ASV.ID"],
  asv.taxa.df[which(asv.taxa.df$ASV.ID %in% colnames(asv.tab.filt) & 
                      asv.taxa.df$Order == "Chloroplast"),"ASV.ID"],
  asv.taxa.df[which(asv.taxa.df$ASV.ID %in% colnames(asv.tab.filt) & 
                      asv.taxa.df$Kingdom == "Eukaryota"),"ASV.ID"]
))

asv.tab.filt2 <- asv.tab.filt[,which(!colnames(asv.tab.filt) %in% excluded.seqs)]
asv.tab.filt2 <- asv.tab.filt2[which(rowSums(asv.tab.filt2) > 0),]

# Supp figure: final read counts post taxa filtering
ggplot(meta[which(!meta$Lab.control),], aes(Sample.type, Reads.postprocessing, color = Cestode.index))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Cestode.index), size = 3, position = position_jitterdodge())+
  scale_x_discrete(labels = c("gut", "wash", "body"))+
  labs(x = "Sample type", y = "Reads post processing", color = "Cestode index")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
        legend.position = "right")
if(!file.exists("figures/Supp_additional_fig_readsPOSTprocessing.png")){
  ggsave("figures/Supp_additional_fig_readsPOSTprocessing.png", units = "in", dpi = 600, width = 6, height = 4)
}


#### final asv table for publication ####
asv.final <- as.data.frame(t(asv.tab))
asv.final$ASV.ID <- row.names(asv.final)
asv.final <- merge(asv.final, asv.taxa.df[,c(1:9)], by = "ASV.ID")
asv.final$ASV_status <- sapply(asv.final$ASV.ID, function(x){
  if(x %in% pcr.false.pos){
    "excluded_PCR_false_positive"
  } else if(x %in% unique(c(contams.blanks, contams.freq))){
    "excluded_contaminant"
  } else if(x %in% c(taxa.mockseqs,"seq179")){
    "mock_community"
  } else if(x %in% colnames(asv.tab.filt2)){
    "included"
  } else if(x %in% excluded.seqs){
    "excluded_eukaryotic_origin"
  } else {
    NA
  }
})
table(asv.final$ASV_status)
asv.final <- merge(asv.final, asv.seqs, by = "ASV.ID")

# Supp table S3
if(!file.exists("figures/Supp_additional_table_ASV_counts_and_taxonomy.csv")){
  # only write file if it doesn't exist
  write.csv(asv.final, "figures/Supp_additional_table_ASV_counts_and_taxonomy.csv",
            row.names = FALSE, quote = FALSE)
}

#### Phyloseq of filtered data ####
## Phyloseq of filtered data ##
ps.filt <- phyloseq(otu_table(asv.tab.filt2, taxa_are_rows=FALSE),
                    tax_table(asv.taxa[colnames(asv.tab.filt2),]),
                    sample_data(meta[which(meta$SampleID.R %in% rownames(asv.tab.filt2)),])
)

#### QC: blanks post filtering ####
ps.blanks <- subset_samples(ps.filt, Lab.neg)
ps.blanks <- filter_taxa(ps.blanks, function(x) sum(x) > 0, TRUE)
ps.blanks
# 3 taxa, 4 samples
otu_table(ps.blanks)
tax_table(ps.blanks)
# blank 2 has Mycoplasma (seq2) and Photbacterium phosphoreum (seq6)
# ddH2O has Photobacterium (seq20)
# most likely cross contamination, not concerned.
plot_bar(ps.blanks, fill = "Genus")+
  theme_classic()
# can exclude blanks from rest of analysis

#### QC: PCR replicates & merging ####
## Ordination ##
ps.samples <- subset_samples(ps.filt, !Lab.control)
ps.samples <- filter_taxa(ps.samples, function(x) sum(x) > 0, TRUE)
ps.samples
# 116 taxa, 134 samples
ps.relab <- transform_sample_counts(ps.samples, function(OTU) OTU/sum(OTU))

ps.clr <- transform_sample_counts(ps.samples, function(OTU) OTU+1) #pseudo-count
ps.clr <- microbiome::transform(ps.clr, 'clr')

if(file.exists("data/ordination_pcr_replicates_clr_nmds.rds")){
  ord.reps <- readRDS("data/ordination_pcr_replicates_clr_nmds.rds")
} else {
  ord.reps <- ordinate(ps.clr, method = "NMDS", distance = "euclidean", k = 3)
  saveRDS(ord.reps, "data/ordination_pcr_replicates_clr_nmds.rds")
}

ord.reps$stress
# 0.12072

# Supp figure: NMDS of PCR replicates
plot_ordination(ps.clr, ord.reps, type="samples", axes=c(1,2), color="Sample.type", shape="PCR.Replicate")+
  geom_vline(xintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_hline(yintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_line(aes(group = SampleID), size = 0.8)+
  geom_point(size=2)+
  scale_x_continuous(limits = c(-20,20))+
  scale_y_continuous(limits = c(-20,20))+
  scale_color_manual(values = gg_color_hue(3), labels = c("salmon gut","cestode wash","cestode body"))+
  coord_fixed(ratio = 1)+
  labs(color="Sample type", shape="PCR replicate")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "right",
        axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
        panel.border = element_rect(size = 1.2))
if(!file.exists("figures/Supp_additional_fig_NMDS_PCR_replicates.png")){
  ggsave("figures/Supp_additional_fig_NMDS_PCR_replicates.png", units = "in", dpi = 600, width = 5, height = 5)
}

# fairly clear grouping by Sample.type type, with cestode wash between salmon gut & cestode body
# replicates also seem to be grouping together - can combine

## merge replicates ##
ps.merge <- merge_samples(ps.samples, "SampleID", )
# warnings generated but seemed to work

# sample_data didn't seem to merge properly
ps.merge <- phyloseq(otu_table(ps.merge),
                     tax_table(ps.merge),
                     sample_data(meta[which(meta$PCR.Replicate == "R1" & !meta$Lab.control),c(2:3,7:14)]))
sample_data(ps.merge)$Reads.sequenced.sum <- sapply(sample_data(ps.merge)$SampleID, function(x){
  sum(meta[which(meta$SampleID == x),"Reads.sequenced"])
})

psm.relab <- transform_sample_counts(ps.merge, function(OTU) OTU/sum(OTU))
psm.clr <- transform_sample_counts(ps.merge, function(OTU) OTU+1) #pseudo-count
psm.clr <- microbiome::transform(psm.clr, 'clr')

# dataframe version for when phyloseq isn't useful
df.merge <- as.data.frame(otu_table(ps.merge))
meta.merged <- data.frame(sample_data(psm.relab))

dfm.relab <- as.data.frame(prop.table(as.matrix(df.merge), margin = 1))
dfm.relab$SampleID <- row.names(dfm.relab)
dfm.relab <- merge(dfm.relab, meta.merged, by = "SampleID")
dfm.relab.m <- melt(dfm.relab, id.vars = grep("seq",names(dfm.relab),invert = T), measure.vars = grep("seq",names(dfm.relab)),
                    value.name = "proportion", variable.name = "ASV.ID")

# Supp figure: rarefaction/accumulation curve to check whether sequencing effort was sufficient
rcurve.all <- rcurve(ps.merge, subsamp = 10^c(1:5), trim = TRUE, add_sample_data = TRUE)
ggplot(rcurve.all, aes(Reads, Richness, group = SampleID, color = Sample.type))+
  geom_point(size = 2)+
  geom_line()+
  scale_x_continuous(trans = "log10")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
        legend.position = "none")+
  facet_wrap(~Sample.type, nrow = 3,
             labeller=labeller(Sample.type = c(salmon_gut="Salmon gut", cestode_wash="Cestode wash", cestode_body="Cestode body")))
if(!file.exists("figures/Supp_additional_fig_accumulation_curve.png")){
  ggsave("figures/Supp_additional_fig_accumulation_curve.png", units = "in", dpi = 600, width = 6, height = 6)
}

#### Ordinations ####
## Sample type
if(file.exists("data/ordination_sample_type_clr_nmds.rds")){
  ord.samples <- readRDS("data/ordination_sample_type_clr_nmds.rds")
} else {
  ord.samples <- ordinate(psm.clr, method = "NMDS", distance = "euclidean", k = 3)
  saveRDS(ord.samples, "data/ordination_sample_type_clr_nmds.rds")
}

ord.samples$stress
"0.1169566"

sample_data(psm.clr)$Sample.type <- factor(sample_data(psm.clr)$Sample.type, levels = c("salmon_gut","cestode_wash","cestode_body"))
p.ord3.1 <- plot_ordination(psm.clr, ord.samples, type="samples", axes=c(1,2), color="Sample.type")+
  geom_point(size=1)+
  scale_x_continuous(limits = c(-20,20))+
  scale_y_continuous(limits = c(-20,20))+
  coord_fixed(ratio = 1)+
  theme_bw()+ 
  theme(text = element_text(color = "black", size = 7), panel.grid = element_blank(), legend.position = "none",
        axis.text = element_text(size = 6, colour = "black"), 
        panel.border = element_rect(size = 0.5))
p.ord3.2 <- plot_ordination(psm.clr, ord.samples, type="samples", axes=c(1,3), color="Sample.type")+
  geom_point(size=1)+
  scale_x_continuous(limits = c(-20,20))+
  scale_y_continuous(limits = c(-20,20))+
  coord_fixed(ratio = 1)+
  theme_bw()+ 
  theme(text = element_text(color = "black", size = 7), panel.grid = element_blank(), legend.position = "none",
        axis.text = element_text(size = 6, colour = "black"), 
        panel.border = element_rect(size = 0.5))
p.ord3.l <- get_legend(
  plot_ordination(psm.clr, ord.samples, type="samples", axes=c(1,2), color="Sample.type")+
    geom_point(size=1)+
    scale_color_manual(values = gg_color_hue(3), labels = c("Salmon gut","Cestode wash","Cestode body"))+
    theme_bw()+ 
    labs(color="Sample type")+
    theme(text = element_text(color = "black", size = 7), panel.grid = element_blank(), legend.position = "right",
          axis.text = element_text(size = 6, colour = "black"), 
          panel.border = element_rect(size = 0.5))
)

# Figure 1 E-F: NMDS of sample type
if(!file.exists("figures/Fig1ef_NMDS_by_sample_type.pdf")){
  pdf("figures/Fig1ef_NMDS_by_sample_type.pdf", 
      width = 6.5, height = 4, colormodel = "srgb", paper = "a4r")
  plot_grid(p.ord3.1, p.ord3.2, p.ord3.l, nrow = 1, rel_widths = c(1,1,0.5))
  dev.off()
}


# PERMANOVA
ord.samples.dist <- phyloseq::distance(psm.clr, method = "euclidean")
adonis(ord.samples.dist ~ Sample.type, data = meta.merged) ###################### 14/2/22 getting weird error?
"            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sample.type  2    4653.6 2326.80  8.7097 0.21395  0.001 ***
Residuals   64   17097.6  267.15         0.78605           
Total       66   21751.2                 1.00000"
betadisper(ord.samples.dist, meta.merged$Sample.type)
"Average distance to median:
  salmon_gut cestode_wash cestode_body 
       16.75        16.25        12.68"
permutest(betadisper(ord.samples.dist, meta.merged$Sample.type))
"          Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
Groups     2  215.07 107.533 6.5706    999  0.003 **
Residuals 64 1047.41  16.366
# difference may be due to difference in dispersions, particularly for cestode body"

## Salmon gut
psm.gut <- subset_samples(ps.merge, Sample.type=="salmon_gut")
psm.gut <- filter_taxa(psm.gut, function(x) sum(x) > 0, TRUE)

psm.gut.relab <- transform_sample_counts(psm.gut, function(OTU) OTU/sum(OTU))
psm.gut.clr <- transform_sample_counts(psm.gut, function(OTU) OTU+1) #pseudo-count
psm.gut.clr <- microbiome::transform(psm.gut.clr, 'clr')

if(file.exists("data/ordination_salmon_gut_clr_nmds.rds")){
  ord.gut <- readRDS("data/ordination_salmon_gut_clr_nmds.rds")
} else {
  ord.gut <- ordinate(psm.gut.clr, method = "NMDS", distance = "euclidean", k = 2)
  saveRDS(ord.gut, "data/ordination_salmon_gut_clr_nmds.rds")
}

ord.gut$stress
"0.1504631"

# Supp figure: NMDS of salmon gut by cestode index
if(!file.exists("figures/Supp_FigS3_NMDS_salmon_gut_by_cestode.pdf")){
  pdf("figures/Supp_FigS3_NMDS_salmon_gut_by_cestode.pdf", width = 3.4, height = 3.4)
  plot_ordination(psm.gut.clr, ord.gut, type="samples", axes=c(1,2), color="Cestode.index")+
    geom_point(size=1)+
    labs(color = "Cestode index")+
    scale_x_continuous(limits = c(-25,25))+
    scale_y_continuous(limits = c(-25,25))+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(text = element_text(color = "black", size = 7), panel.grid = element_blank(), 
          legend.position = "right",
          axis.text = element_text(size = 6, colour = "black"), 
          panel.border = element_rect(size = 0.5))
  dev.off()
}


#### Figure 1: microbiota composition between sample types ####
## top 12 genera (+ dominant Mycoplasmataceae ASV)
psm.genus <- rank_summary(ps.merge, "Genus")
psm.genus$SampleID <- row.names(psm.genus)
psm.genus <- psm.genus[,c(ncol(psm.genus),1:(ncol(psm.genus)-2))]
names(psm.genus)
asv.taxa.df[which(asv.taxa.df$ASV.ID == "seq9"),c("Family","Genus")]
all.equal(psm.genus$Unknown, df.merge$seq9) #TRUE
names(psm.genus)[6] <- "Mycoplasmataceae_sp"
psm.genus$Unassigned <- rowSums(psm.genus[,grep("Unknown",names(psm.genus))])
psm.genus <- psm.genus[,grep("Unknown", names(psm.genus), invert = T)]
psm.genus.top <- top_taxa(psm.genus, sample_var = "SampleID", top = 12)
psm.genus.top <- melt(psm.genus.top, id.vars = "SampleID", value.name = "count", variable.name = "taxon")
psm.genus.top <- merge(psm.genus.top, meta.merged, by = "SampleID", all.x = T, all.y = F)
psm.genus.top$Fish.ID.of <- factor(psm.genus.top$Fish.ID, levels = unique(meta.merged[order(meta.merged$Cestode.index),"Fish.ID"]), ordered = TRUE)
psm.genus.top$Sample.type <- factor(psm.genus.top$Sample.type, levels = c("salmon_gut","cestode_wash","cestode_body"))

# Figure 1B-D: 16S genus composition
if(!file.exists("figures/Fig1bcd_genus_composition_by_sample_type.pdf")){
  pdf("figures/Fig1bcd_genus_composition_by_sample_type.pdf",
      width = 6.5, height = 5, colormodel = "srgb", paper = "a4r")
  ggplot(psm.genus.top, aes(Fish.ID.of, count, fill=taxon))+
    facet_wrap(Sample.type~., nrow = 3, 
               labeller=labeller(Sample.type = c(salmon_gut="Salmon gut", cestode_wash="Cestode wash", cestode_body="Cestode body"))) +
    geom_col(position="fill", width = 1)+
    labs(x="Fish ID", y="Proportion of 16S reads", fill="Genus")+
    scale_fill_manual(values=c(brewer.pal(12, "Set3"), "black","#969696"))+
    geom_vline(xintercept = 9.5, linetype=2)+
    geom_vline(xintercept = 18.5, linetype=2)+
    geom_vline(xintercept = 27.5, linetype=2)+
    theme(text = element_text(color = "black", size = 7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 90),
          legend.position = "right", strip.text = element_text(face = "bold"))
  dev.off()
}


## shared and unique ASVs between sample types (using ps_venn from MicEco package)
## use plot = FALSE to return a list of shared and unique taxa
# presence/absence of taxa, no abundance filter
ps_venn(psm.relab, "Sample.type", weight = FALSE)
ps_venn(psm.gut, "Cestode.present", weight = FALSE)

sample_data(psm.relab)$SampleType.CestodeBin <- sapply(sample_data(psm.relab)$SampleID, function(x){
  df <- data.frame(sample_data(psm.relab))
  if(df[which(df$SampleID == x),"Sample.type"] == "salmon_gut"){
    paste0(c("gut","cestode",df[which(df$SampleID == x),"Cestode.present"]), collapse = "_")
  } else {
    as.character(df[which(df$SampleID == x),"Sample.type"])
  }
})

# colours: gut = "#F8766D", wash = "#00BA38", body = "#619CFF"; cestode = "#fcdf03", no cestode = "#530087"

if(!file.exists("figures/Supp_FigS4c_Venn_16Staxa.png")){
  png("figures/Supp_FigS4c_Venn_16Staxa.png", units = "in", res = 600, width = 6, height = 4)
  ps_venn(psm.relab, "SampleType.CestodeBin", weight = FALSE, fill = NA, col = c("#619CFF","#00BA38","#530087","#fcdf03"), lwd = 2,
          labels = list(labels=NA, font=1), legend = TRUE)
  dev.off()
}


#### Maaslin2 analysis to identify significantly changing genera and ASVs ####

maas.asv.relab <- as.data.frame(otu_table(psm.relab))
maas.gen.relab <- as.data.frame(prop.table(as.matrix(psm.genus[,2:ncol(psm.genus)]), margin = 1))
maas.gen.relab <- maas.gen.relab[,names(maas.gen.relab) != "Unassigned"]

meta.merged$Size.class <- factor(meta.merged$Size.class, levels = c("Small","Medium","Large"), ordered = T)
meta.merged$Sample.type <- factor(meta.merged$Sample.type, levels = c("salmon_gut","cestode_wash","cestode_body"))
meta.merged$Sampling.Date <- as.factor(meta.merged$Sampling.Date)
meta.merged$Samples.per.fish <- sapply(as.character(meta.merged$Fish.ID), function(x) { nrow(meta.merged[which(meta.merged$Fish.ID == x),]) })

# don't re-run if out files already exist

## Genus relab
# Sample.type comparison
if(file.exists("maaslin2_out/GEN_relab_v1_sample_type/all_results.tsv")){
  maas.gen.relab.v1 <- read.delim("maaslin2_out/GEN_relab_v1_sample_type/all_results.tsv")
} else {
  maas.gen.relab.v1 <- Maaslin2(
    input_data = maas.gen.relab[which(meta.merged$Samples.per.fish == 3),],
    input_metadata = meta.merged[which(meta.merged$Samples.per.fish == 3),],
    fixed_effects = c("Sample.type"),
    random_effects = c("Fish.ID","Sampling.Date"),
    reference = c("Sample.type,salmon_gut"),
    output = "maaslin2_out/GEN_relab_v1_sample_type",
    normalization = "NONE", transform = "NONE"
  )
  maas.gen.relab.v1 <- maas.gen.relab.v1$results
}

# subset to salmon gut
if(file.exists("maaslin2_out/GEN_relab_v2_gut_cestode_present/all_results.tsv")){
  maas.gen.relab.v2a <- read.delim("maaslin2_out/GEN_relab_v2_gut_cestode_present/all_results.tsv")
} else {
  maas.gen.relab.v2a <- Maaslin2(
    input_data = maas.gen.relab[which(meta.merged$Sample.type == "salmon_gut"),],
    input_metadata = meta.merged[which(meta.merged$Sample.type == "salmon_gut"),],
    fixed_effects = c("Cestode.present","Feed.Type","Sex","Size.class"),
    random_effects = c("Sampling.Date"),
    reference = c("Size.class,Small"),
    output = "maaslin2_out/GEN_relab_v2_gut_cestode_present",
    normalization = "NONE", transform = "NONE"
  )
  maas.gen.relab.v2a <- maas.gen.relab.v2a$results
}

if(file.exists("maaslin2_out/GEN_relab_v2_gut_cestode_index/all_results.tsv")){
  maas.gen.relab.v2b <- read.delim("maaslin2_out/GEN_relab_v2_gut_cestode_index/all_results.tsv")
} else {
  maas.gen.relab.v2b <- Maaslin2(
    input_data = maas.gen.relab[which(meta.merged$Sample.type == "salmon_gut"),],
    input_metadata = meta.merged[which(meta.merged$Sample.type == "salmon_gut"),],
    fixed_effects = c("Cestode.index","Feed.Type","Sex","Size.class"),
    random_effects = c("Sampling.Date"),
    reference = c("Size.class,Small"),
    output = "maaslin2_out/GEN_relab_v2_gut_cestode_index",
    normalization = "NONE", transform = "NONE"
  )
  maas.gen.relab.v2b <- maas.gen.relab.v2b$results
}

## ASV relab
# Sample.type
if(file.exists("maaslin2_out/ASV_relab_v1_sample_type/all_results.tsv")){
  maas.asv.relab.v1 <- read.delim("maaslin2_out/ASV_relab_v1_sample_type/all_results.tsv")
} else {
  maas.asv.relab.v1 <- Maaslin2(
    input_data = maas.asv.relab[which(meta.merged$Samples.per.fish == 3),],
    input_metadata = meta.merged[which(meta.merged$Samples.per.fish == 3),],
    fixed_effects = c("Sample.type"),
    random_effects = c("Fish.ID","Sampling.Date"),
    reference = c("Sample.type,salmon_gut"),
    output = "maaslin2_out/ASV_relab_v1_sample_type",
    normalization = "NONE", transform = "NONE"
  )
  maas.asv.relab.v1 <- maas.asv.relab.v1$results
}

# subset to fish gut
if(file.exists("maaslin2_out/ASV_relab_v2_gut_cestode_present/all_results.tsv")){
  maas.asv.relab.v2a <- read.delim("maaslin2_out/ASV_relab_v2_gut_cestode_present/all_results.tsv")
} else {
  maas.asv.relab.v2a <- Maaslin2(
    input_data = maas.asv.relab[which(meta.merged$Sample.type == "salmon_gut"),],
    input_metadata = meta.merged[which(meta.merged$Sample.type == "salmon_gut"),],
    fixed_effects = c("Cestode.present","Feed.Type","Sex","Size.class"),
    random_effects = c("Sampling.Date"),
    reference = c("Size.class,Small"),
    output = "maaslin2_out/ASV_relab_v2_gut_cestode_present",
    normalization = "NONE", transform = "NONE"
  )
  maas.asv.relab.v2a <- maas.asv.relab.v2a$results
}

if(file.exists("maaslin2_out/ASV_relab_v2_gut_cestode_index/all_results.tsv")){
  maas.asv.relab.v2b <- read.delim("maaslin2_out/ASV_relab_v2_gut_cestode_index/all_results.tsv")
} else {
  maas.asv.relab.v2b <- Maaslin2(
    input_data = maas.asv.relab[which(meta.merged$Sample.type == "salmon_gut"),],
    input_metadata = meta.merged[which(meta.merged$Sample.type == "salmon_gut"),],
    fixed_effects = c("Cestode.index","Feed.Type","Sex","Size.class"),
    random_effects = c("Sampling.Date"),
    reference = c("Size.class,Small"),
    output = "maaslin2_out/ASV_relab_v2_gut_cestode_index",
    normalization = "NONE", transform = "NONE"
  )
  maas.asv.relab.v2b <- maas.asv.relab.v2b$results
}

# Supp Table: Taxa associations
maas.gen.relab.v1[which(maas.gen.relab.v1$qval < 0.2),]
maas.gen.relab.v2a[which(maas.gen.relab.v2a$qval < 0.2),]
#none
maas.gen.relab.v2b[which(maas.gen.relab.v2b$qval < 0.2),]

maas.asv.relab.v1[which(maas.asv.relab.v1$qval < 0.2),]
maas.asv.relab.v1$ASV_genus <- sapply(maas.asv.relab.v1$feature, function(x){
  if(is.na(asv.taxa.df[as.character(x),"Genus"])){
    asv.taxa.df[as.character(x),"Family"]
  } else {
    asv.taxa.df[as.character(x),"Genus"]
  }
})

maas.asv.relab.v2a[which(maas.asv.relab.v2a$qval < 0.2),]
#none
maas.asv.relab.v2b[which(maas.asv.relab.v2b$qval < 0.25),]
maas.asv.relab.v2b$ASV_genus <- sapply(maas.asv.relab.v2b$feature, function(x){
  if(is.na(asv.taxa.df[as.character(x),"Genus"])){
    asv.taxa.df[as.character(x),"Family"]
  } else {
    asv.taxa.df[as.character(x),"Genus"]
  }
})

maas.asv.relab.sample_type.sigseqs <- unique(
  maas.asv.relab.v1[which(maas.asv.relab.v1$qval < 0.2),"feature"]
)

maas.asv.relab.cestode.sigseqs <- unique(
  maas.asv.relab.v2b[which(maas.asv.relab.v2b$qval < 0.25),"feature"]
)

maas.asv.relab.df <- subset(
  dfm.relab.m, ASV.ID %in% c(maas.asv.relab.sample_type.sigseqs,maas.asv.relab.cestode.sigseqs)
)

maas.asv.relab.df$Genus <- sapply(maas.asv.relab.df$ASV.ID, function(x){
  if(is.na(asv.taxa.df[as.character(x),"Genus"])){
    asv.taxa.df[as.character(x),"Family"]
  } else {
    asv.taxa.df[as.character(x),"Genus"]
  }
})
maas.asv.relab.df$ASV.ID.of <- factor(
  maas.asv.relab.df$ASV.ID, levels = unique(maas.asv.relab.df[order(maas.asv.relab.df$Genus),"ASV.ID"]), ordered = TRUE
)

# Supp table of Maaslin2 results
if(!file.exists("figures/Supp_TableS1_maaslin2_results.csv")){
  maas.gen.relab.v1$model <- "Genus_m1"
  maas.gen.relab.v2b$model <- "Genus_m2"
  maas.gen.relab.v1$ASV_genus <- NA
  maas.gen.relab.v2b$ASV_genus <- NA
  
  maas.asv.relab.v1$model <- "ASV_m1"
  maas.asv.relab.v2b$model <- "ASV_m2"
  
  maas.res <- rbind.data.frame(
    maas.gen.relab.v2b[which(maas.gen.relab.v2b$qval < 0.25),c(2,1,3:9,11,10)],
    maas.asv.relab.v2b[which(maas.asv.relab.v2b$qval < 0.25),c(2,1,3:11)],
    maas.gen.relab.v1[which(maas.gen.relab.v1$qval < 0.25),c(2,1,3:9,11,10)],
    maas.asv.relab.v1[which(maas.asv.relab.v1$qval < 0.25),c(2,1,3:11)]
  )
  write.csv(maas.res, "figures/Supp_TableS1_maaslin2_results.csv", row.names = F, quote = F)
}


#### Top ASV and genera plots ####

# Supp Figure: Cestode index (ASV)
p.bx.asv <-
ggplot(maas.asv.relab.df[which(maas.asv.relab.df$ASV.ID %in% maas.asv.relab.cestode.sigseqs &
                                 maas.asv.relab.df$Sample.type == "salmon_gut"),],
       aes(Cestode.index, proportion))+
  facet_wrap(~Genus + ASV.ID, scales = "free_y", nrow = 1)+
  geom_boxplot(outlier.colour = NA, fill = NA)+
  geom_point(aes(color=Cestode.index), size = 2, position = position_jitterdodge())+
  labs(x="Cestode index", y="Proportion of 16S reads")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))


# Supp Figure: Cestode index/present (genus)
psm.genus.top.rel <- top_taxa(psm.genus, sample_var = "SampleID", top = 12)
psm.genus.top.rel <- as.data.frame(prop.table(as.matrix(psm.genus.top.rel[,2:ncol(psm.genus.top.rel)]), margin = 1))
psm.genus.top.rel$SampleID <- row.names(psm.genus.top.rel)
psm.genus.top.rel <- melt(psm.genus.top.rel, id.vars = "SampleID", value.name = "proportion", variable.name = "taxon")
psm.genus.top.rel <- merge(psm.genus.top.rel, meta.merged, by = "SampleID", all.x = T, all.y = F)
psm.genus.top.rel$Sample.type <- factor(psm.genus.top.rel$Sample.type, levels = c("salmon_gut","cestode_wash","cestode_body"))

p.bx.gen.1 <-
  ggplot(psm.genus.top.rel[which(psm.genus.top.rel$Sample.type == "salmon_gut" &
                                   psm.genus.top.rel$taxon %in% c("Mycoplasma","Photobacterium")),],
         aes(Cestode.index, proportion))+
  facet_wrap(~taxon) +
  geom_boxplot(outlier.colour = NA, fill = NA)+
  geom_point(aes(color=Cestode.index), size = 2, position = position_jitterdodge())+
  scale_y_continuous(limits = c(0,1))+
  labs(x="Cestode index", y="Proportion of 16S reads")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

p.bx.gen.2 <-
  ggplot(psm.genus.top.rel[which(psm.genus.top.rel$Sample.type == "salmon_gut" &
                                   psm.genus.top.rel$taxon %in% c("Carnobacterium","Aliivibrio","Brevinema")),],
         aes(Cestode.present, proportion))+
  facet_wrap(~taxon) +
  geom_boxplot(outlier.colour = NA, fill = NA)+
  geom_point(aes(color=Cestode.present), size = 2, position = position_jitterdodge())+
  scale_color_manual(values = c("#440154FF","#12b3b3"))+
  scale_y_continuous(limits = c(0,1))+
  labs(x="Cestode present", y="Proportion of 16S reads")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

# Supp Figure: Cestode index/presence (genus + ASV combined)
p.bx.gen <- 
  plot_grid(p.bx.gen.1, p.bx.gen.2, align = "h", axis = "tblr", 
            nrow = 1, rel_widths = c(1,1.2), labels = c("a","b"))

plot_grid(p.bx.gen, p.bx.asv, align = "hv", axis = "tblr", nrow = 2, labels = c("","c"), rel_heights = c(1,0.8))
if(!file.exists("figures/Supp_FigS5_cestode_comb_genusASV_boxplot.pdf")){
  ggsave("figures/Supp_FigS5_cestode_comb_genusASV_boxplot.pdf", width = 6.5, height = 6.5)
}

# Supp Figure: Sample type (genus + ASV combined)
maas.plot.st.p1 <-
  ggplot(psm.genus.top.rel[which(psm.genus.top.rel$taxon 
                                 %in% c("Photobacterium","Carnobacterium","Aliivibrio",
                                        "Ureaplasma","Mycoplasma","Brevinema",
                                        "Lactobacillus","Mycoplasmataceae_sp")),],
         aes(Sample.type, proportion))+
  facet_wrap(~taxon, nrow = 2)+
  geom_boxplot(outlier.colour = NA, fill = NA)+
  geom_point(aes(color=Sample.type), size = 2, position = position_jitterdodge())+
  scale_x_discrete(labels = c("gut","wash","body"))+
  labs(x="Sample type", y="Proportion of 16S reads")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

maas.plot.st.p2 <- 
  ggplot(maas.asv.relab.df[which(maas.asv.relab.df$ASV.ID %in% maas.asv.relab.sample_type.sigseqs),],
         aes(Sample.type, proportion))+
  facet_wrap(~Genus + ASV.ID.of, scales = "free_y", nrow = 4)+
  geom_boxplot(outlier.colour = NA, fill = NA)+
  geom_point(aes(color=Sample.type), size = 2, position = position_jitterdodge())+
  scale_x_discrete(labels = c("gut","wash","body"))+
  labs(x="Sample type", y="Proportion of 16S reads")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

plot_grid(maas.plot.st.p1, maas.plot.st.p2, align = "hv", axis = "tblr", nrow = 2, rel_heights = c(0.75,1), labels = c("a","b"))
if(!file.exists("figures/Supp_FigS6_sampletype_comb_genusASV_boxplot.pdf")){
  ggsave("figures/Supp_FigS6_sampletype_comb_genusASV_boxplot.pdf", units = "in", dpi = 600, width = 10, height = 10)
}


#### Alpha diversity ####
alpha.asv <- data.frame(SampleID=row.names(df.merge),
                        Richness=specnumber(df.merge),
                        Shannon=diversity(t(df.merge), index = "shannon"))

alpha.asv <- merge(alpha.asv, meta.merged, by = "SampleID")
alpha.asv$Sample.type <- factor(alpha.asv$Sample.type, levels = c("salmon_gut","cestode_wash","cestode_body"))

# Supp Figure: alpha diversity
p.alpha1 <-
ggplot(alpha.asv, aes(Sample.type,Richness))+ #
  facet_grid(~Cestode.index)+
  geom_boxplot(outlier.color = NA)+
  geom_point(aes(color = Sample.type), size = 2, position = position_jitterdodge())+
  scale_x_discrete(labels = c("gut","wash","body"))+
  labs(x = "Sample type", y = "Richness (no. of ASVs)")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

summary(lm(Richness ~ Cestode.index, alpha.asv, subset = Sample.type=="salmon_gut"))
"Residual standard error: 7.733 on 26 degrees of freedom
Multiple R-squared:  0.06923,	Adjusted R-squared:  -0.03817 
F-statistic: 0.6446 on 3 and 26 DF,  p-value: 0.5933"

summary(lm(Richness ~ Cestode.index + Reads.sequenced.sum, alpha.asv, subset = Sample.type=="salmon_gut"))
"Residual standard error: 7.328 on 25 degrees of freedom
Multiple R-squared:  0.1962,	Adjusted R-squared:  0.06758 
F-statistic: 1.525 on 4 and 25 DF,  p-value: 0.2252"

summary(lm(Richness ~ Sample.type, alpha.asv))
"Residual standard error: 5.883 on 64 degrees of freedom
Multiple R-squared:   0.17,	Adjusted R-squared:  0.1441 
F-statistic: 6.555 on 2 and 64 DF,  p-value: 0.002573"

summary(lm(Richness ~ Sample.type + Reads.sequenced.sum, alpha.asv))
"Residual standard error: 5.736 on 63 degrees of freedom
Multiple R-squared:  0.2233,	Adjusted R-squared:  0.1863 
F-statistic: 6.036 on 3 and 63 DF,  p-value: 0.001113"

## Alpha diversity using Hill's numbers framework
alpha.asv.hill <- data.frame(SampleID=row.names(otu_table(psm.relab)),
                             Hill.q1=hill_taxa(otu_table(psm.relab), q = 1),
                             Hill.q2=hill_taxa(otu_table(psm.relab), q = 2))
alpha.asv <- merge(alpha.asv, alpha.asv.hill, by = "SampleID")
p.alpha2 <-
ggplot(alpha.asv, aes(Sample.type,Hill.q1))+
  facet_grid(~Cestode.index)+
  geom_boxplot(outlier.color = NA)+
  geom_point(aes(color = Sample.type), size = 2, position = position_jitterdodge())+
  scale_x_discrete(labels = c("gut","wash","body"))+
  labs(x = "Sample type", y = "Shannon diversity")+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

plot_grid(p.alpha1, p.alpha2, align = "hv", axis = "tblr", nrow = 2, labels = c("a","b"))
if(!file.exists("figures/Supp_FigS4ab_alpha_diversity_both.png")){
  ggsave("figures/Supp_FigS4ab_alpha_diversity_both.png", units = "in", dpi = 600, width = 6.5, height = 5)
}

ggplot(alpha.asv, aes(Sample.type,Hill.q1))+
  geom_boxplot(outlier.color = NA)+
  geom_point(aes(color = Sample.type), size = 2, position = position_jitterdodge())+
  theme_bw()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        legend.position = "none", strip.text = element_text(face = "bold"))

summary(lm(Hill.q1 ~ Cestode.index, alpha.asv, subset = Sample.type=="salmon_gut"))
"Residual standard error: 3.137 on 26 degrees of freedom
Multiple R-squared:  0.08566,	Adjusted R-squared:  -0.01984 
F-statistic: 0.8119 on 3 and 26 DF,  p-value: 0.4988"

summary(lm(Hill.q1 ~ Cestode.index + Reads.sequenced.sum, alpha.asv, subset = Sample.type=="salmon_gut"))
"Residual standard error: 3.121 on 25 degrees of freedom
Multiple R-squared:  0.1298,	Adjusted R-squared:  -0.009477 
F-statistic: 0.9319 on 4 and 25 DF,  p-value: 0.4614"

summary(lm(Hill.q1 ~ Sample.type, alpha.asv))
"Residual standard error: 2.352 on 64 degrees of freedom
Multiple R-squared:  0.06119,	Adjusted R-squared:  0.03185 
F-statistic: 2.086 on 2 and 64 DF,  p-value: 0.1326"
# cestode body significantly lower than gut, whereas cestode wash not

summary(lm(Hill.q1 ~ Sample.type + Reads.sequenced.sum, alpha.asv))
"Residual standard error: 2.364 on 63 degrees of freedom
Multiple R-squared:  0.06638,	Adjusted R-squared:  0.02192 
F-statistic: 1.493 on 3 and 63 DF,  p-value: 0.225"
# ok, tendancy now only "weak" evidence for cestode body to be lower in diversity than gut (p=0.06)

#### Fig 2 & 3: Mycoplasma clades based on 16S tree ####
myco.salmon <- c("seq91","seq2","seq3","seq124","seq64","seq131","seq170","seq817","seq26","seq156")
myco.cestode1 <- c("seq29","seq73","seq18","seq21","seq361","seq1")
myco.cestode2 <- c("seq7","seq10","seq58","seq57")

myco.df <- data.frame(
  SampleID=row.names(df.merge),
  Myco_Salmon=rowSums(df.merge[,myco.salmon]),
  Myco_Cestode1=rowSums(df.merge[,myco.cestode1]),
  Myco_Cestode2=rowSums(df.merge[,myco.cestode2]),
  Total=rowSums(df.merge)
)
myco.df.m <- melt(myco.df, id.vars = c("SampleID","Total"), measure.vars = c(2:4),
                  value.name = "count", variable.name = "Myco_clade")
myco.df.m$relab <- myco.df.m$count / myco.df.m$Total
myco.df.m <- merge(myco.df.m, meta.merged, by = "SampleID", all.x = T, all.y = F)

if(!file.exists("figures/Fig3_mycoplasma_clade_by_sample_type.pdf")){
  pdf("figures/Fig3_mycoplasma_clade_by_sample_type.pdf", 
      width = 6.5, height = 2.5, colormodel = "srgb", paper = "a4r")
  ggplot(myco.df.m[grep("Myco",myco.df.m$Myco_clade),], aes(Sample.type, relab))+
    facet_wrap(~Myco_clade, nrow = 1,
               labeller=labeller(Myco_clade = c(Myco_Salmon="Salmon clade", Myco_Cestode1="Cestode clade 1", Myco_Cestode2="Cestode clade 2")))+
    geom_line(aes(group = Fish.ID, color = Cestode.index, linetype = as.factor(Samples.per.fish)))+
    geom_point(aes(color = Cestode.index), size = 2)+
    scale_x_discrete(labels = c("salmon\ngut","cestode\nwash","cestode\nbody"))+
    scale_linetype_manual(values = c(0,3,1))+
    labs(x="Sample type", y="Proportion of 16S reads", color = "Cestode index", linetype = "Samples per fish")+
    theme_bw()+
    theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
          line = element_line(size = 0.5),
          axis.text = element_text(size = 6, colour = "black"),
          legend.position = "right", strip.text = element_text(face = "bold"))
  dev.off()
}


myco.tree.order <- read.delim("data/mycoplasma_16S_tree_ASV_tip_order.txt", header = F)
myco.tree.order <- grep("seq",myco.tree.order$V1, value = T)
myco.tree.samples <- meta.merged[order(meta.merged$Sample.type, meta.merged$Cestode.index),"SampleID"]

dfm.myco.relab <- dfm.relab.m[which(dfm.relab.m$ASV.ID %in% myco.tree.order),]
dfm.myco.relab[which(dfm.myco.relab$proportion == 0),"proportion"] <- NA
dfm.myco.relab$ASV.ID <- factor(dfm.myco.relab$ASV.ID, levels = rev(myco.tree.order))
dfm.myco.relab$SampleID <- factor(dfm.myco.relab$SampleID, levels = myco.tree.samples)

if(!file.exists("figures/Fig2b_mycoplasma_abundance_presence.pdf")){
  pdf("figures/Fig2b_mycoplasma_abundance_presence.pdf",
      width = 10, height = 6, colormodel = "srgb", paper = "a4r")
  ggplot(dfm.myco.relab, aes(SampleID, ASV.ID))+
    geom_point(aes(color = Sample.type, size = proportion), alpha = 0.75)+
    geom_hline(yintercept = 9.5, color = "grey50", linetype = 2)+
    geom_hline(yintercept = 19.5, color = "grey50", linetype = 2)+
    labs(x = "", y = "")+
    theme_classic()+
    theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
          line = element_line(size = 0.5),
          axis.text = element_text(size = 6, colour = "black"),
          legend.position = "right", axis.text.x = element_blank())
  dev.off()
}


#### Size/weight vs Cestode index ####
## HoloFish entire cohort
# Cestode distribution by feed type
holofish.cestode.distrib <- as.data.frame(prop.table(table(meta.holofish[,c("Feed.Type","Cestode.index")]), margin = 1))
holofish.cestode.distrib$y_pos <- NA
holofish.cestode.distrib[holofish.cestode.distrib$Feed.Type == "Feed1","y_pos"] <- 
  cumsum(rev(holofish.cestode.distrib[holofish.cestode.distrib$Feed.Type == "Feed1","Freq"]))
holofish.cestode.distrib[holofish.cestode.distrib$Feed.Type == "Feed2","y_pos"] <- 
  cumsum(rev(holofish.cestode.distrib[holofish.cestode.distrib$Feed.Type == "Feed2","Freq"]))
holofish.cestode.distrib$Freq.2 <- round(holofish.cestode.distrib$Freq, digits = 3)

ggplot(holofish.cestode.distrib, aes(Feed.Type, Freq, fill = Cestode.index))+
  geom_col(position = "stack")+
  scale_fill_manual(values = viridis::viridis(4))+
  labs(x = "Feed type", y = "Proportion of individuals", fill = "Cestode index")+
  theme_classic()+
  theme(text = element_text(color = "black", size = 7), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5), axis.text = element_text(size = 6, colour = "black"),
        legend.position = "right")
if(!file.exists("figures/Supp_FigS1_cestode_index_vs_feedtype.png")){
  ggsave("figures/Supp_FigS1_cestode_index_vs_feedtype.png", units = "in", dpi = 600,  width = 3, height = 3)
}


# stats
table(meta.holofish$Cestode.present)
"FALSE  TRUE 
   85   378 "
prop.table(table(meta.holofish$Cestode.present))
"    FALSE      TRUE 
0.1835853 0.8164147"

prop.table(table(meta.holofish$Cestode.index))
"        0          1          2          3 
0.18358531 0.41468683 0.31533477 0.08639309 "

prop.table(table(meta.holofish[,c("Feed.Type","Cestode.index")]), margin = 1)

summary(lm(Gutted.Weight.kg ~ Cestode.index + Feed.Type + Sex, meta.holofish))
"Residual standard error: 1.213 on 456 degrees of freedom
  (1 observation deleted due to missingness)
Multiple R-squared:  0.1761,	Adjusted R-squared:  0.1671 
F-statistic:  19.5 on 5 and 456 DF,  p-value: < 2.2e-16"

## Cestode study subset
meta.ind <- subset(meta.holofish, Fish.ID %in% unique(meta$Fish.ID))
table(meta.ind[,c("Size.class","Cestode.present")])
"          Cestode.present
Size.class FALSE TRUE
    Large      2    5
    Medium     6    6
    Small      1   10"

prop.table(table(meta.ind[,c("Size.class","Cestode.present")]), margin = 2)
"         Cestode.present
Size.class     FALSE      TRUE
    Large  0.2222222 0.2380952
    Medium 0.6666667 0.2857143
    Small  0.1111111 0.4761905"
table(meta.ind[,c("Feed.Type","Cestode.present")])
"         Cestode.present
Feed.Type FALSE TRUE
    Feed1     5   10
    Feed2     4   11"
table(meta.ind[,c("Sex","Cestode.present")])
"   Cestode.present
Sex FALSE TRUE
  F     3    5
  M     6   16"
prop.table(table(meta.ind[,c("Sex","Cestode.present")]), margin = 1)
"   Cestode.present
Sex     FALSE      TRUE
  F 0.3750000 0.6250000
  M 0.2727273 0.7272727"

summary(lm(Gutted.Weight.kg ~ Cestode.index + Feed.Type + Sex, meta.ind))
"Residual standard error: 1.177 on 24 degrees of freedom
Multiple R-squared:  0.2147,	Adjusted R-squared:  0.05109 
F-statistic: 1.312 on 5 and 24 DF,  p-value: 0.2921"

## Supp Figure: index vs weight ##
sfig.kg.a <- ggplot(meta.holofish, aes(Cestode.index, Gutted.Weight.kg))+
  geom_violin(draw_quantiles = 0.5)+
  geom_point(aes(color = Cestode.index), size = 2, position = position_jitterdodge())+
  labs(x = "Cestode index", y = "Gutted weight (kg)", title = "HoloFish cohort")+
  theme_classic()+
  theme(text = element_text(color = "black", size = 11), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5), axis.text = element_text(size = 12, colour = "black"),
        legend.position = "none")

sfig.kg.b <- ggplot(meta.ind, aes(Cestode.index, Gutted.Weight.kg))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Cestode.index), size = 2, position = position_jitterdodge())+
  labs(x = "Cestode index", y = "Gutted weight (kg)", title = "Cestode investigation")+
  theme_classic()+
  theme(text = element_text(color = "black", size = 11), rect = element_rect(size = 0.5),
        line = element_line(size = 0.5), axis.text = element_text(size = 12, colour = "black"),
        legend.position = "none")

if(!file.exists("figures/Supp_FigS2_cestode_index_vs_weight.pdf")){
  pdf("figures/Supp_FigS2_cestode_index_vs_weight.pdf", width = 8, height = 4)
  plot_grid(sfig.kg.a, sfig.kg.b, align = "hv", axis = "tblr", labels = c("a","b"), nrow = 1)
  dev.off()
}

### MAG functional comparison ####
gene.clusters <- read.delim("data/mycoplasma_MAG_gene_clusters.txt")

gene.clusters$annotated <- sapply(seq(1,nrow(gene.clusters),1), function(i){
  k <- gene.clusters[i,"KOfam"] == ""
  c <- gene.clusters[i,"COG20_FUNCTION_ACC"] == ""
  p <- gene.clusters[i,"Pfam"] == ""
  if(k & c & p){
    #all blank
    FALSE
  } else {
    TRUE
  }
})

gc.id.lists <- list(
  Ml = unique(gene.clusters[which(gene.clusters$genome_name == "Candidatus_Mycoplasma_lavaretus"),"gene_cluster_id"]),
  Msm = unique(gene.clusters[which(gene.clusters$genome_name == "Candidatus_Mycoplasma_salmoninae_mykiss"),"gene_cluster_id"]),
  Mss = unique(gene.clusters[which(gene.clusters$genome_name == "Candidatus_Mycoplasma_salmoninae_salar"),"gene_cluster_id"]),
  M_iowae = unique(gene.clusters[which(gene.clusters$genome_name == "Mycoplasma_iowae"),"gene_cluster_id"]),
  M_penetrans = unique(gene.clusters[which(gene.clusters$genome_name == "Mycoplasma_penetrans"),"gene_cluster_id"]),
  M_mobile = unique(gene.clusters[which(gene.clusters$genome_name == "Mycoplasma_mobile"),"gene_cluster_id"]),
  CE_seq1 = unique(gene.clusters[which(gene.clusters$genome_name == "CE_seq1"),"gene_cluster_id"]),
  CE_seq7 = unique(gene.clusters[which(gene.clusters$genome_name == "CE_seq7"),"gene_cluster_id"]),
  CE_seq10 = unique(gene.clusters[which(gene.clusters$genome_name == "CE_seq10"),"gene_cluster_id"]),
  U_diversum = unique(gene.clusters[which(gene.clusters$genome_name == "Ureaplasma_diversum"),"gene_cluster_id"]),
  U_urealyticum = unique(gene.clusters[which(gene.clusters$genome_name == "Ureaplasma_urealyticum"),"gene_cluster_id"])
)


### diff versions
ggVennDiagram(gc.id.lists[c("CE_seq7","CE_seq1","Mss","Msm","Ml")], 
              category.names = c("CE_seq7","CE_seq1","Ms. salar","Ms. mykiss", "M. lavaretus"),
              label = "count", label_alpha = 0)+
  scale_fill_gradient(low = NA, high = NA)+
  scale_color_manual(values = c("#BC93FF","#A4EEDE","#FFCDAB","#FFE6D5","#E5D5FF"))+
  theme(legend.position = "none")

ggVennDiagram(gc.id.lists[c("CE_seq7","CE_seq1","Mss","Msm","Ml")], 
              category.names = c("CE_seq7","CE_seq1","Ms. salar","Ms. mykiss", "M. lavaretus"),
              label = "count", label_alpha = 0)+
  scale_fill_gradient(low = NA, high = NA)+
  scale_color_manual(values = rep("black",5))+
  theme(legend.position = "none")
if(!file.exists("figures/Supp_FigS7a_Venn_MAG_geneclusters.png")){
  ggsave("figures/Supp_FigS7a_Venn_MAG_geneclusters.png", units = "in", dpi = 600, width = 6.5, height = 6)
}

table(gene.clusters[,c("genome_name","annotated")])
prop.table(table(gene.clusters[,c("genome_name","annotated")]), margin = 1)


#### END ####






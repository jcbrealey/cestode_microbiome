#!/usr/bin/Rscript

### Rscript for processing silva reference sequences for Mycoplasma tree
## v2 23/10/2020 - cutting down on reference sequences

setwd("/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/Qiime2_phylogeny")

# Mycoplasma
ltp <- read.delim("../SILVA_LTP_132_SSU/LTPs132_SSU_mycoplasma_IDs.txt", header = F)
rep <- read.delim("Mycoplasma_representative_strains_201009.txt", sep=",", header = T)
silva <- read.delim("silva_nr_v138_custom_201008_Mycoplasma_species_IDs.txt", sep=",", header = F)

names(silva) <- c("ID.full","Genus","Strain")
silva$ID.short <- sapply(as.character(silva$ID.full), function(x){
  strsplit(x,"\\.")[[1]][1]
})
silva$Species <- sapply(as.character(silva$Strain), function(x){
  strsplit(x," ")[[1]][1]
})
silva$LTP.in <- silva$ID.short %in% ltp$V1

myco.sp <- data.frame(Species=unique(silva$Species))

myco.sp$ID.use <- sapply(as.character(myco.sp$Species), function(x){
  df <- silva[which(silva$Species == x),]
  if(nrow(df) == 1) {
    # if only 1 entry, use that ID
    as.character(df$ID.full)
  } else {
    if(nrow(df[which(df$LTP.in),]) > 0){
      # if in LTP, use first ID
      as.character(df[which(df$ID.short %in% intersect(df$ID.short, ltp$V1))[1],"ID.full"])
    } else {
      # pull rep strain from list and pick first instance 
      as.character(df[grep(rep[which(rep$species == x),"Strain"], df$Strain)[1],"ID.full"])
    }
  }
})

# exclude contaminated moatsii sequence
myco.sp <- myco.sp[grep("moatsii",myco.sp$Species,invert=T),]
# exclude Limborg sequences as will add in separately
myco.sp <- myco.sp[grep("Limborg",myco.sp$ID.use,invert=T),]
# exclude LQ020882 "Mycoplasma pneumonia" as an error in taxonomy
myco.sp <- myco.sp[which(myco.sp$ID.use != "LQ020882.8.1460"),]

# check for key species 
print(paste("iowae:", grep("iowae", myco.sp$Species, value = T)))
print(paste("penetrans:", grep("penetrans", myco.sp$Species, value = T)))
print(paste("mobile:", grep("mobile", myco.sp$Species, value = T)))
print(paste("No. of Mycoplasma refs:", nrow(myco.sp)))

## only include subset of reference sequences
jacob.refs.myco <- c("collis","molare","conjunctivae","bovoculi","ovipneumoniae","hyopneumoniae",
                "flocculare","dispar","hyorhinis","pulmonis","crocodyli","synoviae","glycophilum",
                "canis","cynos","columbinum","meleagridis","lipofaciens","fermentans",
                "agalactiae","bovis","mobile","phocidae","arginini","hominis","arthritidis",
                "cloacale","mycoides","capricolum","leachii","putrefacians","yeatsii",
                "gallisepticum","genitalium","pneumoniae","alvi","pirum","iowae","penetrans",
                "girerdii","haemofelis","suis","wenyonii","insons","cavipharyngis")
myco.sp <- myco.sp[which(myco.sp$Species %in% jacob.refs.myco),]
print(paste("No. of Mycoplasma refs v2:", nrow(myco.sp)))

# Ureaplasma
ltp <- read.delim("../SILVA_LTP_132_SSU/LTPs132_SSU_ureaplasma_IDs.txt", header = F)
silva <- read.delim("silva_nr_v138_custom_201008_Ureaplasma_species_IDs.txt", sep=",", header = F)
# no rep strains needed
names(silva) <- c("ID.full","Genus","Strain")
silva$ID.short <- sapply(as.character(silva$ID.full), function(x){
  strsplit(x,"\\.")[[1]][1]
})
silva$Species <- sapply(as.character(silva$Strain), function(x){
  strsplit(x," ")[[1]][1]
})
silva$LTP.in <- silva$ID.short %in% ltp$V1

urea.sp <- data.frame(Species=unique(silva$Species))

urea.sp$ID.use <- sapply(as.character(urea.sp$Species), function(x){
  df <- silva[which(silva$Species == x),]
  if(nrow(df) == 1) {
    # if only 1 entry, use that ID
    as.character(df$ID.full)
  } else {
    if(nrow(df[which(df$LTP.in),]) > 0){
      # if in LTP, use first ID
      as.character(df[which(df$ID.short %in% intersect(df$ID.short, ltp$V1))[1],"ID.full"])
    } else {
      # miroungigenitalium should be only one, pick first instance
      as.character(df[1,"ID.full"])
    }
  }
})

print(paste("No. of Ureaplasma refs:", nrow(urea.sp)))

## only include subset of reference sequences
jacob.refs.urea <- c("diversum","urealyticum")
urea.sp <- urea.sp[which(urea.sp$Species %in% jacob.refs.urea),]
print(paste("No. of Ureaplasma refs v2:", nrow(urea.sp)))

# Save IDs
writeLines(c(myco.sp$ID.use, urea.sp$ID.use), con = "silva_nr_v138_custom_201023_TREE_IDS.txt")



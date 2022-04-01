#!/bin/bash -l

# Pull out ref sequences for Mycoplasma family from silva

REF_file=reference_databases/silva_nr_v138_custom_200821_species_assignment.fa.gz
WDDIR=/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/Qiime2_phylogeny
OUT_base=silva_nr_v138_custom_201008

# mycoplasma
zcat $REF_file | grep " Mycoplasma " --no-group-separator -A1 > $WDDIR/${OUT_base}_Mycoplasma_species.fa

grep "^>" $WDDIR/${OUT_base}_Mycoplasma_species.fa > $WDDIR/${OUT_base}_Mycoplasma_species_IDs.txt
sed -i 's/>//g' $WDDIR/${OUT_base}_Mycoplasma_species_IDs.txt
sed -i 's/ Mycoplasma /,Mycoplasma,/g' $WDDIR/${OUT_base}_Mycoplasma_species_IDs.txt

# ureaplasma
zcat $REF_file | grep " Ureaplasma " --no-group-separator -A1 > $WDDIR/${OUT_base}_Ureaplasma_species.fa

grep "^>" $WDDIR/${OUT_base}_Ureaplasma_species.fa > $WDDIR/${OUT_base}_Ureaplasma_species_IDs.txt
sed -i 's/>//g' $WDDIR/${OUT_base}_Ureaplasma_species_IDs.txt
sed -i 's/ Ureaplasma /,Ureaplasma,/g' $WDDIR/${OUT_base}_Ureaplasma_species_IDs.txt

# Pull out ref IDs from silva Living Tree Project
REFDIR2=/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/SILVA_LTP_132_SSU
grep "Mycoplasma " $REFDIR2/LTPs132_SSU_compressed.fasta | cut -f1 > $REFDIR2/LTPs132_SSU_mycoplasma_IDs.txt
sed -i 's/>//g' $REFDIR2/LTPs132_SSU_mycoplasma_IDs.txt

grep "Ureaplasma " $REFDIR2/LTPs132_SSU_compressed.fasta | cut -f1 > $REFDIR2/LTPs132_SSU_ureaplasma_IDs.txt
sed -i 's/>//g' $REFDIR2/LTPs132_SSU_ureaplasma_IDs.txt



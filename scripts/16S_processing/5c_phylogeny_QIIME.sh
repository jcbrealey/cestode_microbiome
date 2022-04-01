#!/bin/bash -l

# Phylogeny of input sequences using latest version of QIIME2 (2020)
# alignment with mafft
# mask positions with gaps in >20% of sequences (--p-max-gap-frequency 0.2)
# tree with raxml and bootstrapping


source activate qiime2-2020.8

# Set up data
DATADIR=/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/Qiime2_phylogeny
ff=mycoplasmaFAM_ASV_refs_201023.fna

cd $DATADIR

# merge reference sequences and ASVs if haven't already
#module load seqtk
#seqtk subseq silva_nr_v138_custom_201008_Mycoplasma_species.fa silva_nr_v138_custom_201023_TREE_IDS.txt > mycoplasma_refs_201023.fna
#seqtk subseq silva_nr_v138_custom_201008_Ureaplasma_species.fa silva_nr_v138_custom_201023_TREE_IDS.txt > ureaplasma_refs_201023.fna

#cat mycoplasmaFAM_ASV_201008.fna mycoplasma_ref_extras_GenBank.fna mycoplasma_refs_201023.fna ureaplasma_refs_201023.fna > $ff
#zcat ../SILVA_138_SSU_dada2/Limborg_16S_seqs_species.fa.gz | grep "Mycoplasma" -A1 --no-group-separator >> $ff
#zcat ../SILVA_138_SSU_dada2/silva_species_assignment_v138.fa.gz | grep ">AY876289" -A1 --no-group-separator >> $ff

# Run QIIME commands

qiime tools import --input-path $ff --output-path ${ff%fna}qza --type 'FeatureData[Sequence]'
qiime alignment mafft --i-sequences ${ff%fna}qza --o-alignment ${ff%.fna}_aln.qza
qiime alignment mask --i-alignment ${ff%.fna}_aln.qza --p-max-gap-frequency 0.2 --p-min-conservation 0.4 --o-masked-alignment ${ff%.fna}_aln_masked.qza
qiime phylogeny raxml-rapid-bootstrap --i-alignment ${ff%.fna}_aln.qza --p-seed 248 --p-rapid-bootstrap-seed 625 --p-bootstrap-replicates 1000 --p-substitution-model GTRGAMMA --p-n-threads 5 --o-tree ${ff%.fna}_aln_masked_tree.qza
qiime tools export --input-path ${ff%.fna}_aln_masked_tree.qza --output-path .
mv tree.nwk ${ff%.fna}_aln_masked.tree

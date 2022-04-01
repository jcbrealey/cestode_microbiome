Scripts for processing 16S data, from raw (demultiplexed) Illumina sequences to ASVs with DADA2.

Raw sequences can be downloaded from ENA.

Note scripts were run on a local server and inputs/file structures are set up accordingly (i.e. they don't follow the Git repository structure).

1. adapter removal: 1_adapter_removal.sh
2. DADA2 filtering: 2_DADA2_filtering.R
3. DADA2 infer ASV: 3_DADA2_inferASVs.R
4. assign taxonomy: 4_DADA2_assignTaxonomy.R
5. phylogeny:
	5a. Prepare refs1: 5a_phylogeny_get_ref_seqIDs.sh
	5b. Prepare refs2: 5b_phylogeny_process_ref_seqs.R
	5c. Run phylogeny: 5c_phylogeny_QIIME.sh

Reference databases required for input:
reference_databases/silva_nr_v138_custom_200821_train_set.fa.gz"
reference_databases/silva_nr_v138_custom_200821_species_assignment.fa.gz"
reference_databases/LTPs132_SSU_compressed.fasta
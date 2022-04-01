### In progress (scripts haven't been added yet) ###
Scripts for processing shotgun metagenomics data, from raw (demultiplexed) BGI sequences to MAGs with Anvi'o.

Raw sequences can be downloaded from ENA.

Note scripts were run on a local server and inputs/file structures are set up accordingly (i.e. they don't follow the Git repository structure).

QC and Processing:
1. adapter removal: adrm_pe_210302_repeat_test.sh
2. remove PCR duplicates: rmdup_210303_repeat.sh
3. re-match PE reads: repair_210303.sh
4. filter out phiX reads: rmphix_210308.sh
5. filter out salmon reads: rmsalmon_210311.sh
6. filter out cestode reads: rmtapeworm_210313.sh
7. filter out human reads: rmhuman_210314.sh

Reference genomes for processing:
phiX: NC_001422.1.phiX174.fna 
salmon: GCF_000233375.1_ICSASG_v2_genomic.fna
cestode: Cestode_genomes_201110/bwa_cestode_genomes_n24/cestode_genomes_n24_201118.fna
human: GCF_000001405.39_GRCh38.p13_genomic.fna

MAG assembly:
6. coassembly of samples: megahit_coassembly_210315.sh
7. map samples back to contigs: contig_mapping_210316.sh
8. profile with anvi'o: mag_recovery_step2-4_210316_COMG_minC1000.sh

MAG phylogeny
## TBD


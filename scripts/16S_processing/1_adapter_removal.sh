#!/bin/bash -l

## Preprocessing of 16S tapeworm reads
## DADA2/Begum approach
## Step 1: AdapterRemoval 
##         - Keep PE reads separate
##         - remove Illumina adapters + primers from 3' ends
##         - going to do quality trimming and length filtering in DADA2
##         - generate useful stats


module load AdapterRemoval
echo $(module list)

# adapters - should match the END of fwd and rev reads, respectively
# Not including 16S primer region (will do in DADA2)
ad1="GTAGTCCCTGTCTCTTATACACATCTCCGAGCCCACGAGAC........ATCTCGTATGCCGTCTTCTGCTTG"
ad2="CTGTCTCTTATACACATCTGACGCTGCCGACGA........GTGTAGATCTCGGTGGTCGCCGTATCATT"

DATADIR=RAW_DATA
OUTDIR=/groups/hologenomics/jaelleb/TAPEWORM_16S_200729/P1_adrm_200811_PE

mkdir $OUTDIR
cd $DATADIR

for i in *R1_001.fastq.gz;
do
    echo $i
    cd $OUTDIR
    sname=$(echo $i | sed -e 's/TOG-2JTMM-\(.*\)_L.*/\1/')
    if [[ "$sname" == "Mock"* || "$name" == "Blank"* || "$name" == "ddH2O"* ]]; 
    then
       bname=${sname}_adrm
    else
       bname=i${sname}_adrm #add i to fish sample nr so that files do not start with a number
    fi
    AdapterRemoval --file1 $DATADIR/$i --file2 $DATADIR/${i%1_001.fastq.gz}2_001.fastq.gz --basename $bname --adapter1 $ad1 --adapter2 $ad2 --gzip --threads 4;
done

# generate stats
cd $OUTDIR
echo "#sample,preprocessed,adapter-containing,PE.truc,singletons,ave.length.bp,fwd.16S.match,rev.16S.match" >> log_adrm_stats_200811.txt

for i in *settings;
do
    sid=$(echo "${i%_adrm.settings}")
    c1=$(grep "Total number of read pairs:" $i | cut -f2 --delimiter ':')
    c2=$(grep "Number of reads with adapters" $i | cut -f2 --delimiter ':')
    let c3=$(zcat ${i%settings}pair1.truncated.gz | wc -l)/4
    let c4=$(zcat ${i%settings}singleton.truncated.gz | wc -l)/4
    l1=$(grep "Average length of retained reads:" $i | cut -f2 --delimiter ':')
    c5=$(zcat ${i%settings}pair1.truncated.gz | grep -c -i ^ccta.ggg..gca.cag)
    c6=$(zcat ${i%settings}pair2.truncated.gz | grep -c -i ^ggactac..gggtatctaat)
    echo "$sid,$c1,$c2,$c3,$c4,$l1,$c5,$c6" >> log_adrm_stats_200811.txt;
done

# remove unnecessary files
rm *.settings
rm *.discarded.gz
rm *.singleton.truncated.gz

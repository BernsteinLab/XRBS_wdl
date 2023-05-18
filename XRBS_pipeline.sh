#!/bin/bash -l
#$ -N XRBS_pipeline
#$ -cwd
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=8G
#$ -e XRBS_pipeline.log
#$ -o XRBS_pipeline.log
#$ -j yes
#$ -l h_rt=36:00:00


'''
Pipeline to:
- Trim adapters from fastq and quality trim
- Convert fastq to a 3-letter genome
- Align to 3-letter genome
- Convert alignment
- Mark duplicates, sort and index
- Call methylated bases

Have to change the path manually for each new run
The script has two inputs: the sample name and the date of the run,
which is assigned manually to separate the different
analyses

date and folder structures have to be created
manually, outside of this script

Best way to run this analysis:
Create the date directory with the fastq, trim and align subdirs
Prepare a list of the samples to be analyzed
and run the script in a for loop as in:

while read sample; do
  qsub /seq/epiprod02/sgaldon/Alba/scripts/XRBS_pipeline.sh $sample path
done <sample_names.txt
'''


source /broad/software/scripts/useuse
reuse -q .root-v5.34.01
reuse -q UGER
conda activate python2.7
use Samtools
use .openssl-1.0.2g
#use .tabix-0.2.6
use Java-1.8
use .picard-tools
use .bwa-0.7.4
use .fastqc-0.11.9
use .r-3.6.0-bioconductor

sample=$1
path=$2

fastq_path="$path"/fastqs/
trim_path="$path"/trim/
align_path="$path"/align/

# echo $sample

mkdir "$fastq_path"/fastqc 

fastqc "$fastq_path"/"$sample"_R1.fastq.gz "$fastq_path"/"$sample"_R2.fastq.gz -o "$fastq_path"/fastqc 

## Run cutadapt

~/.local/bin/cutadapt --discard -a GCTCTTCCGATCT -o "$trim_path"/"$sample"_R2.cutadapt.fastq.gz -p "$trim_path"/"$sample"_R1.cutadapt.fastq.gz "$fastq_path"/"$sample"_R2.fastq.gz "$fastq_path"/"$sample"_R1.fastq.gz > "$trim_path"/"$sample"_R1.cutadapt_report.txt

/seq/epiprod02/hovestadt/software/TrimGalore-0.6.5/trim_galore --path_to_cutadapt ~/.local/bin/cutadapt -o "$trim_path" --paired --illumina --nextseq 20 "$trim_path"/"$sample"_R1.cutadapt.fastq.gz "$trim_path"/"$sample"_R2.cutadapt.fastq.gz

#Run Methyltools
#https://github.com/hovestadt/methylCtools
#fqconv script converts illumina BS sequencing reads to 3-letter alphabet.

/seq/epiprod02/hovestadt/software/methylctools-0.9.4/methylCtools fqconv -1 "$trim_path"/"$sample"_R1.cutadapt_val_1.fq.gz -2 "$trim_path"/"$sample"_R2.cutadapt_val_2.fq.gz "$align_path"/"$sample".conv.fq

bwa mem -t 4 -p -M -T 0 /seq/epiprod02/hovestadt/genomes/methylCtools-GRCh38/GRCh38.conv.fa "$align_path"/"$sample".conv.fq | samtools view -Sb - > "$align_path"/"$sample".reads.conv.bam


#Run MethylCtools bconv
#reads are converted back to original state. positions are stored in read id.
/seq/epiprod02/hovestadt/software/methylctools-0.9.4/methylCtools bconv "$align_path"/"$sample".reads.conv.bam -m "$align_path"/"$sample".human.conv.sort.metrics.txt - | samtools sort -T "$align_path"/"$sample".human.sort -@ 4 - > "$align_path"/"$sample".human.sort.bam

samtools index "$align_path"/"$sample".human.sort.bam

# java -Xmx4g -jar $PICARD MarkDuplicates I="$align_path"/"$sample".human.sort.bam O="$align_path"/"$sample".human.sort.dup.bam M="$align_path"/"$sample".human.sort.dup.bam.metrics TMP_DIR='/seq/epiprod02/sgaldon/tmp'
# samtools index "$align_path"/"$sample".human.sort.dup.bam &
# samtools flagstat "$align_path"/"$sample".human.sort.dup.bam > "$align_path"/"$sample".human.sort.dup.bam.flagstat &
# wait

samtools flagstat "$align_path"/"$sample".human.sort.bam > "$align_path"/"$sample".human.sort.bam.flagstat

Rscript /seq/epiprod02/sgaldon/Alba/scripts/filter.vh20200112.R "$align_path"/"$sample".human.sort.bam

/seq/epiprod02/hovestadt/software/methylctools-0.9.4/methylCtools bcall --trimPE --metrics "$align_path"/"$sample".human.sort.filter.call.metrics /seq/epiprod02/hovestadt/genomes/methylCtools-GRCh38/GRCh38.pos.gz "$align_path"/"$sample".human.sort.bam.filter.bam - | bgzip > "$align_path"/"$sample".human.sort.filter.call.gz

/seq/epiprod02/hovestadt/software/methylctools-0.9.4/methylCtools bcall --trimPE --metrics "$align_path"/"$sample".human.sort.filter.CH.call.metrics /seq/epiprod02/hovestadt/genomes/methylCtools-GRCh38/GRCh38.CHpos_chr1.pos.gz "$align_path"/"$sample".human.sort.bam.filter.bam - | bgzip > "$align_path"/"$sample".human.sort.filter.CH.call.gz

samtools flagstat "$align_path"/"$sample".human.sort.bam.filter.bam > "$align_path"/"$sample".human.sort.bam.filter.bam.flagstat


gzip "$align_path"/"$sample".conv.fq

##############
## Calculate coverage

# samtools depth "$align_path"/"$sample".human.sort.dup.subsampled.bam | gzip > "$align_path"/"$sample".human.sort.dup.depth.tsv.gz

# awk '{ key=($1 > $2 ? $1 FS $2 : $2 FS $1) }
#     NR==FNR { a[key]; next }
#     key in a
# ' <(zcat "$align_path"/"$sample".human.sort.dup.subsampled.call.gz) <(zcat "$align_path"/"$sample".human.sort.dup.depth.tsv.gz) > "$align_path"/"$sample".human.sort.dup.subsampled.subsetDepth.tsv.gz

# rm "$align_path"/"$sample".human.sort.dup.depth.tsv.gz

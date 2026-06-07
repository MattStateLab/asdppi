#!/bin/bash

#SBATCH --job-name=FOXP1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kelsey.hennick@ucsf.edu
#SBATCH --cpus-per-task=16
#SBATCH --mem=50gb
#SBATCH --time=20:00:00
#SBATCH --output=/nowakowskilab/data1/kelsey_new/kelsey/20260214_FOXP1-R513H_CUTTag/cat/logs/FOXP1_%j.log
#SBATCH -p cpu_medium


# Global Functions ---------------------------------

#Remove and remake file
rmMk()
{
rm -rf $1
mkdir $1
}


# Global variables ------------------------------------------------

module load CBI miniconda3
conda activate cuttagenv2

# variables
srcFileList="/nowakowskilab/data1/kelsey_new/kelsey/20260214_FOXP1-R513H_CUTTag/cat/fileLists/FOXP1_fileList.txt"
srcFileLine=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $srcFileList)
srcFile=$(echo ${srcFileLine} | awk '{print $1}')
ctrlFile=$(echo ${srcFileLine} | awk '{print $2}')

# destination folders 
destParent="/nowakowskilab/data1/kelsey_new/kelsey/20260214_FOXP1-R513H_CUTTag/cat/outputs/FOXP1"
mkdir=${destParent}/${srcFile}
dest="${destParent}/${srcFile}"

# make folders
mkdir ${dest}/files
mkdir ${dest}/trim
mkdir ${dest}/bowtie
mkdir ${dest}/picard
mkdir ${dest}/MACS2

destFiles="/nowakowskilab/data1/kelsey_new/kelsey/20260214_FOXP1-R513H_CUTTag/cat/FASTQ"
destTrim="${dest}/trim"
destBowtie="${dest}/bowtie"
destPicard="${dest}/picard"
destMACS="${dest}/MACS2"
#destStarch="${dest}/starch"



# Notebook Chunks ------------------------------------------------------------------------

# Run trimgalore
runTrimGalore(){
rmMk $destTrim
cd $destTrim

module load CBI trimgalore

trim_galore -q 20 --stringency 5 --paired --fastqc --length 20 -o $destTrim ${destFiles}/${srcFile}_R1_001.fastq.gz ${destFiles}/${srcFile}_R2_001.fastq.gz 
trim_galore -q 20 --stringency 5 --paired --fastqc --length 20 -o $destTrim ${destFiles}/${ctrlFile}_R1_001.fastq.gz ${destFiles}/${ctrlFile}_R2_001.fastq.gz

}

# Align to hg38 using Bowtie2
runBowtie(){
rmMk $destBowtie
cd $destBowtie

module load CBI bowtie2

btIndHs="/nowakowskilab/data1/kelsey_new/kelsey/hg38_2/index"

# Logical alignent parameters: I have this page from Harvard Informatics describing bowtie2 for ATACseq alignment settings.
# temprorarily commenting this out
bowtie2 --end-to-end --very-sensitive --phred33 -I 10 -X 700 -x ${btIndHs} -1 ${destTrim}/${srcFile}_R1_001_val_1.fq.gz -2 ${destTrim}/${srcFile}_R2_001_val_2.fq.gz -S ${destBowtie}/${srcFile}.sam &>${destBowtie}/${srcFile}_alignSummary.txt
bowtie2 --end-to-end --very-sensitive --phred33 -I 10 -X 700 -x ${btIndHs} -1 ${destTrim}/${ctrlFile}_R1_001_val_1.fq.gz -2 ${destTrim}/${ctrlFile}_R2_001_val_2.fq.gz -S ${destBowtie}/${ctrlFile}.sam &>${destBowtie}/${ctrlFile}_alignSummary.txt

module purge
}

# Use Picard to sort the alignment and remove duplicates
runPicard(){
rmMk $destPicard
cd $destPicard

#module load CBI picard

picard="java -jar /c4/home/kmhennick/.conda/envs/cuttagenv2/share/picard-3.4.0-0/picard.jar"

$picard SortSam I=${destBowtie}/${srcFile}.sam O=${destPicard}/${srcFile}.sorted.sam SORT_ORDER=coordinate
#a couple extra steps because there is no header in these files for some reason
#SortSame works without, just not MarkDuplicates
module load CBI samtools
samtools view -H ${destPicard}/${srcFile}.sorted.sam | grep '@RG' #to check if there is an RG, if not run below line
samtools addreplacerg -r "@RG\tID:ReadGroup1\tSM:"${srcFile}"\tPL:Illumina\tLB:Library.fa" -o ${destPicard}/${srcFile}.sorted.rg.sam ${destPicard}/${srcFile}.sorted.sam
#markDuplicates is standard
$picard MarkDuplicates I=${destPicard}/${srcFile}.sorted.rg.sam O=${destPicard}/${srcFile}.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=${destPicard}/${srcFile}.rmDup.txt

$picard SortSam I=${destBowtie}/${ctrlFile}.sam O=${destPicard}/${ctrlFile}.sorted.sam SORT_ORDER=coordinate
#added addreplacerg
samtools addreplacerg -r "@RG\tID:ReadGroup1\tSM:"${ctrlFile}"\tPL:Illumina\tLB:Library.fa" -o ${destPicard}/${ctrlFile}.sorted.rg.sam ${destPicard}/${ctrlFile}.sorted.sam
$picard MarkDuplicates I=${destPicard}/${ctrlFile}.sorted.rg.sam O=${destPicard}/${ctrlFile}.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=${destPicard}/${ctrlFile}.rmDup.txt

module purge
}

# MACS 
runMacs(){
rmMk $destMACS
cd $destMACS
mkdir lambda
mkdir nolambda

#cp ${srcDir}/${srcFile}.bam $destMACS

macs2=""
module load CBI samtools
samtools view -S -b ${destPicard}/${srcFile}.rmDup.sam >${destMACS}/${srcFile}.rmDup.bam
samtools view -S -b ${destPicard}/${ctrlFile}.rmDup.sam >${destMACS}/${ctrlFile}.rmDup.bam


macs2 callpeak -t ${destMACS}/${srcFile}.rmDup.bam -c ${destMACS}/${ctrlFile}.rmDup.bam -f BAMPE -g hs -q 0.01 -B --SPMR --keep-dup all --nomodel --nolambda -n $srcFile --outdir ${destMACS}/nolambda

module purge
}


# run chunks -------------------------------------------------------------------------

# Comment out the sections you don't need
runTrimGalore
runBowtie
runPicard
runMacs








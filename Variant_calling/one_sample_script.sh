# Activate the creadted environment (variant_calling)
conda activate variant_calling

## Workflow for single sample:
######################################

## Download the working raw data:
# Create a new folder named raw_data to store the example data in it 
mkdir raw_data && cd raw_data
# Download the example data in the rawdata drirectory
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz 
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz 


## Check the reads Quality control:
fastqc *.fastq.gz
# Move the results output files into a new directory named untrimmed_qc_results 
mkdir -p /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/untrimmed_qc_results 
mv *_fastqc.zip *_fastqc.html /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/untrimmed_qc_results
# After seeing the QC reports, it seems that forward strand has a problem in the sequance quality and sequence content; while only sequence content only in reverse strand
# So that we will do trimming step to enhance the quality control for both
# Run trimming
trimmomatic PE -phred33 ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz\
 ERR8774458_1.trim.fastq.gz ERR8774458_2.trim.fastq.gz\
 ERR8774458_1.unpaired.fastq.gz ERR8774458_2.unpaired.fastq.gz \
  SLIDINGWINDOW:5:30 MINLEN:50 HEADCROP:5
# Check the output trimmed files which should smaller than the input before trimming
ls -l -h ERR8774458*
# Move the trimmed output files into a new dirsctory named trimmed_results 
mkdir -p /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results 
mv *.trim* *.unpaired* /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results 
# Recheck the quality control after the trimming (this part need to be arranged again)
fastqc /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results/*.trim.fastq.gz


## Alignment to Genome:
# Download the reference in a new directory named ref_seq
mkdir /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq 
wget -P /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq https://zenodo.org/records/10886725/files/Reference.fasta 
# Indexing the reference genome
bwa index /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq/Reference.fasta
# Before perform read mapping to the genome, we will set the working variables
ref=/Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq/Reference.fasta
read1=/Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/raw_data/ERR8774458_1.fastq.gz
read2=/Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/raw_data/ERR8774458_2.fastq.gz
# Aligning the reads to the refernce genome 
bwa mem $ref $read1 $read2 > output.sam
# Converting sam to BAM format 
samtools view -S -b output.sam > output.bam
# Sorting BAM file 
samtools sort -o output.sorted.bam output.bam
# Learning more about the statistics of the sorted bam file
samtools flagstat output.sorted.bam
# Assess the alignment (Visualization)
# To visualize the alignment files, we need to index the BAM file
samtools index output.sorted.bam
# Viewing our mapped reads with reference genome
samtools tview output.sorted.bam $ref


## Variant calling:
# Calculate the read coverage of position in the genome
bcftools mpileup -O b -o output.bcf \
-f /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq/Reference.fasta output.sorted.bam
# Detecting the single nucleotice variants (SNVs)
bcftools call --ploidy 1 -m -v -o variants.vcf output.bcf
# Filter and report the SNV cariant in VCF format 
vcfutils.pl varFilter variants.vcf > final_variant.vcf
# Explore the vcf file 
less -S final_variant.vcf
# Explore number of variants
grep -v "#" final_variant.vcf | wc -l

# Activate the creadted environment (variant_calling)
conda activate variant_calling

## Workflow for multiple samples:
#################################

## Download the working raw data:
# Create a new folder named raw_data to store the example data in it 
mkdir raw_datasets && cd raw_datasets
# Download the example data in the rawdata drirectory
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz


## Check the reads Quality control:
fastqc *.fastq.gz
# Move the results output files into a new directory named untrimmed_qc
mkdir -p /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/untrimmed_qc
mv *_fastqc.zip *_fastqc.html /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/untrimmed_qc


# Run trimming
#For loop to trimming all the paired dataset
# For each R1 in the datasets
for R1_file in *_R1.fastq.gz; do
    # Extracting R2 file name 
    R2_file="${R1_file/_R1/_R2}"
    # Perform trimming 
    trimmomatic PE -phred33 "$R1_file" "$R2_file" \
    "${R1_file%.fastq.gz}_trim.fastq.gz" "${R1_file%.fastq.gz}_unpaired.fastq.gz" \
    "${R2_file%.fastq.gz}_trim.fastq.gz" "${R2_file%.fastq.gz}_unpaired.fastq.gz" \
    SLIDINGWINDOW:4:15 MINLEN:36 
done

# Check the output trimmed files which should smaller than the input before trimming
ls -l -h 
# Move the trimmed output files into a new dirsctory named trimmed_results 
mkdir -p /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results 
mv *_trim* *_unpaired* /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results 
# Recheck the quality control after the trimming (this part need to be arranged again)
fastqc /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results/*_trim.fastq.gz


## Alignment to Genome:
# Download the reference in a new directory named ref_seq_all
mkdir /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq_all 
wget -P /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq_all https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta 
# Indexing the reference genome
bwa index /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq_all/Reference.fasta

# Create a new directory  named mapping_files to get sam output files after the mapping 
mkdir -p /Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/mapping_files 
# Before perform read mapping to the genome, we will set the working variables
ref_all=/Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/ref_seq_all/Reference.fasta
input_dir="/Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/trimmed_results"
output_dir="/Users/hadeeribrahiem/Desktop/Bioinfo_online_courses/BioHack_Internship/Stage_2/mapping_files"

# For loop to aligning the reads to the refernce genome 
for R1_file in "$input_dir"/*_R1_trim.fastq.gz; do
    R2_file="${R1_file/_R1/_R2}" 
    #Extract the file name without directory (Get the base filename)
    base="${R1_file##*/}"  
    # Remove the _R1 substring from the filename
    base="${base%_R1_trim.fastq.gz}" 
    # Construct full output file paths
    sam_output="${output_dir}/${base}.sam"
    bam_output="${output_dir}/${base}.bam"
    sorted_bam_output="${output_dir}/${base}.sorted.bam"
    # Mapping
    bwa mem "$ref_all" "$R1_file" "$R2_file" > "$sam_output"
    # Converting sam to BAM format 
    samtools view -S -b "$sam_output" > "$bam_output"
     # Sorting BAM file 
    samtools sort -o "$sorted_bam_output" "$bam_output"
done

# Learning more about the statistics of the sorted bam file
for each_file in "$output_dir"/*.sorted.bam; do
    echo "Statistics for file: "${each_file##*/}" "
    samtools flagstat $each_file
done 

# Assess the alignment (Visualization)
# To visualize the alignment files, we need to index the BAM file
for each_file in "$output_dir"/*.sorted.bam; do 
    samtools index $each_file
    # Viewing our mapped reads with reference genome
    samtools tview $each_file $ref_all
done     


## Variant calling:
# Loop through each sorted BAM file in the directory
for bam_file in "$output_dir"/*.sorted.bam; do
    # Output BCF file
    output_bcf="${bam_file%.sorted.bam}.bcf"
    # Output VCF file
    output_vcf="${bam_file%.sorted.bam}.vcf"
    # Calculate the read coverage of position in the genome
    bcftools mpileup -O b -o "$output_bcf" -f "$ref_all" "$bam_file"
    # Detecting the single nucleotice variants (SNVs)
    bcftools call --ploidy 1 -m -v -o "$output_vcf" "$output_bcf"
   # Filter and report the SNV cariant in VCF format
    vcfutils.pl varFilter "$output_vcf" > "${output_vcf%.vcf}_filtered.vcf"
done


# Explore the vcf file 
for vcf_file in "$output_dir"/*_filtered.vcf; do
    less -S $vcf_file
    # Explore number of variants
    num_variants=$(grep -v "#" $vcf_file | wc -l)
    #Extract the file name without directory (Get the base filename)
    base="${vcf_file##*/}"  
    # Remove the _R1 substring from the filename
    base="${base%_filtered.vcf}" 
    echo "Number of Variants in $base: $num_variants"
done    


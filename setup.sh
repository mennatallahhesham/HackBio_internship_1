## Installing the working tools:
# Add two channels in conda
conda --add channels bioconda
conda --add channels forge
# Confirm existance of the two channels 
conda config --show channels
# Create an environment to install the working tools inside it 
conda create -n variant_calling fastqc multiqc bwa samtools bcftools trimmomatic vcftools

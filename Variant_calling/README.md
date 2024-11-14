# Variant-Calling-Workflow

This repository provides a comprehensive workflow for variant calling analysis using bioinformatics tools. It encompasses both single sample and multiple samples' variant calling processes.

## Usage Instructions

### Environment Setup

Run the provided setup script `setup.sh` to install the necessary tools using Conda.

```bash
bash setup.sh
```
## Work flow 

1- Download the exapmle data (parird-end FASTQ Reads).

2- Download the reference Genome.

3- Assess Quality of Reads.

4- Read Mapping.

5- Filter Alignment Records.

6- Mark Duplicates.

7- Variant calling.

8- Basic statistic.

9- Filter the Variants.

## Tools Used

The following bioinformatics tools are used in this workflow:

- FastQC
- Trimmomatic
- BWA
- SAMtools
- BCFtools
- VCFtools
- MultiQC

## Note

- Make sure to update file paths according to your directory structure.
- Customize the workflow based on your specific requirements.

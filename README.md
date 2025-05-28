# Genepattern-mutect2

# Mutect2 Command Wrapper

This repository provides a small Python-based wrapper to simplify running GATK Mutect2 by handling file indexing and assembling the command-line call. It consists of two main modules:

* **Argument parsing & command builder** (`cmd_builder.py`)
* **Preprocessing helpers** (`preprocess.py`)

## Files

* **preprocess.py**: functions to index FASTA, BAM, and feature files (VCFs).
* **cmd\_builder.py**: constructs the `gatk Mutect2` command based on provided args.
* **run\_mutect2.py**: the main entry point; parses CLI flags, calls preprocessors, then runs Mutect2.

## Usage

### Required arguments

* `-R`, `--reference`: Path to the reference FASTA (e.g. `hg38.fa`)
* `-I`, `--input_bam`: Path to the tumor BAM file
* `-O`, `--output_vcf`: Path to write the output VCF

### Optional arguments

* `-N`, `--normal_bam`: Path to the normal BAM (tumor-normal mode) (Required for somatic varient mode)
* `--tumor_sample`: Tumor sample name label 
* `--normal_sample`: Normal sample name label (Required if -N is passed in)
* `--germline_resource`: Germline resource VCF for allele frequencies
* `--panel_of_normals`: Panel-of-normals VCF
* `-L`, `--intervals`: BED file or intervals string

## Behavior

1. **Index checks**:

   * Reference FASTA (`.fai` & `.dict`)
   * Tumor (and optional normal) BAM (`.bai`)
   * Germline resource & PON VCFs (`.tbi`/`.idx`)
     Existing index files are detected and skipped.

2. **Execution**: Runs GATK Mutect2 with the assembled command.

## Further Documentation

For more details on available Mutect2 parameters and advanced usage, please refer to the official GATK Mutect2 documentation:

* [https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-Mutect2)

# Mutect2 Somatic Variant Calling Module

---

## Overview

This is a GenePattern module wrapping GATK Mutect2. This module supports optional use of a germline reference, panel of normals (PoN), intervals, and parallel scatter–gather execution.

Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations. The caller uses a Bayesian somatic genotyping model that differs from the original MuTect by Cibulskis et al., 2013 and uses the assembly-based machinery of HaplotypeCaller.

---

## Features

- Tumor vs. matched Normal somatic variant calling  
- Support for **Germline Resource** (e.g. gnomAD)  
- Support for **Panel of Normals (PoN)** 
- Optional **interval-based calling** (BED or interval list)  
- Optional **parallel processing** (scatter–gather)  
- Produces **both unfiltered and filtered VCFs**  
- Outputs Mutect2 **filtering and contamination metrics**  

---

## Input Files

### **1. Tumor BAM** *(required)*  
Coordinate-sorted BAM of the tumor sample with `.bai` index.

### **2. Reference Genome FASTA** *(required)*  
Must include `.fai` and `.dict` index files.

### **3. Normal BAM** *(optional but recommended)*  
Matched normal BAM with index.

### **4. Germline Resource (VCF)** *(optional but recommended)*  
Population allele-frequency resource such as `gnomAD`.

### **5. Panel of Normals (PoN)** *(optional but recommended)*  
A PoN VCF identifying recurrent technical artefacts.

### **6. Intervals (BED or .interval_list)** *(optional)*  
Used to restrict calling to target regions.


---

## Output Files

| File | Description |
|------|-------------|
| `<OutputName>.unfiltered.vcf.gz` | Raw candidate variants |
| `<OutputName>.vcf.gz` | Filtered somatic variants |
| `<OutputName>.vcf.gz.tbi` | Tabix index |
| `<OutputName>_metrics.txt` | Filter stats & QC metrics |
| (Optional) intermediate scatter VCFs | Only if parallel mode enabled |

---

## Parameters

| Parameter | Description | Required |
|----------|-------------|----------|
| **Reference** | Reference genome FASTA | Yes |
| **Tumor BAM** | BAM file for tumor sample | Yes |
| **Tumor Sample** | SM tag for tumor | Yes |
| **Output Name** | Prefix for output files | Yes |
| **Normal BAM** | Matched normal BAM | No |
| **Normal Sample** | SM tag for normal | Required if Normal BAM provided |
| **Germline Resource** | gnomAD or similar VCF | No |
| **Panel of Normals** | PoN VCF | No |
| **Intervals** | Regions to call variants | No |
| **Parallel Processing** | "Yes" or "No" | Yes |

---

## Parallel Processing (Scatter–Gather)

If enabled, the module:

1. Splits genomic intervals into chunks  
2. Runs parallel Mutect2 calls  
3. Gathers VCFs into a single unified result

Recommended for WGS or large targeted panels. Can cause boundary defects in output

If **disabled**, Mutect2 runs a single job across all regions.

---

## Sample Data

**Tumor & Normal Reads: FastQ/BAM data for tumor and normal samples from GIAB (HG008, Liss_lab) via FTP:**

https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/

**Reference Genome: hg38 assembly:**

"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

**Germline Resource**

https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz

**Pannel of Normals**

https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz

**Intervals List**

https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list

---

## References

**GATK Mutect2 documentation**
https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

---

## Version

Mutect2 Module — **v1.0**  

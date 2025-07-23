import argparse
import os
from preprocess import has

def parse_args():
    p = argparse.ArgumentParser(description="Run GATK Mutect2 with optional parallel processing.")

    # Required
    p.add_argument("-R", "--reference", required=True, help="Reference genome FASTA")
    p.add_argument("-I", "--input_bam", required=True, help="Tumor BAM file")

    # Optional, but flagged always (default='', nargs='?' to allow --flag without value)
    p.add_argument("-O", "--output_vcf", nargs='?', const='', default='', help="Output VCF filename")
    p.add_argument("-N", "--normal_bam", nargs='?', const='', default='', help="Normal BAM file")
    p.add_argument("--tumor_sample", nargs='?', const='', default='', help="Tumor sample name")
    p.add_argument("--normal_sample", nargs='?', const='', default='', help="Normal sample name")
    p.add_argument("--germline_resource", nargs='?', const='', default='', help="Germline resource VCF")
    p.add_argument("--panel_of_normals", nargs='?', const='', default='', help="Panel of Normals VCF")
    p.add_argument("-L", "--intervals", nargs='?', const='', default='', help="Genomic intervals")
    p.add_argument("--parallel", action="store_true", help="Enable parallel processing across chromosomes")

    # Parse and apply default output if needed
    args = p.parse_args()
    if not args.output_vcf:
        tumor_basename = os.path.basename(args.input_bam)
        tumor_root = os.path.splitext(tumor_basename)[0]
        args.output_vcf = f"{tumor_root}_mutect2.vcf"

    return args



def build_mutect2_cmd(args):
    cmd = ["gatk", "Mutect2",
           "-R", args.reference,
           "-I", args.input_bam]

    if has(args.normal_bam):
        cmd += ["-I", args.normal_bam]
    if has(args.tumor_sample):
        cmd += ["--tumor-sample", args.tumor_sample]
    if has(args.normal_sample):
        cmd += ["--normal-sample", args.normal_sample]
    if has(args.germline_resource):
        cmd += ["--germline-resource", args.germline_resource]
    if has(args.panel_of_normals):
        cmd += ["--panel-of-normals", args.panel_of_normals]
    if hasattr(args, "intervals") and has(args.intervals):
        cmd += ["-L", args.intervals]
    if hasattr(args, "f1r2_tar_gz") and has(args.f1r2_tar_gz):
        cmd += ["--f1r2-tar-gz", args.f1r2_tar_gz]
    if hasattr(args, "bam_output") and has(args.bam_output):
        cmd += ["--bam-output", args.bam_output]
    if hasattr(args, "native_pair_hmm_threads") and args.native_pair_hmm_threads and args.native_pair_hmm_threads > 0:
        cmd += ["--native-pair-hmm-threads", str(args.native_pair_hmm_threads)]

    cmd += ["-O", args.output_vcf]
    return cmd

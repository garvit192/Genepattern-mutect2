import argparse
from preprocess import has

def parse_args():
    p = argparse.ArgumentParser(
        description="Index files & run GATK Mutect2"
    )
    # required arguments
    p.add_argument(
        "-R", "--reference", required=True,
        help="Reference genome FASTA (e.g. hg38.fa)"
    )
    p.add_argument(
        "-I", "--input_bam", required=True,
        help="Input BAM (tumor)"
    )
    p.add_argument(
        "-O", "--output_vcf", required=True,
        help="Output VCF"
    )
    # optional arguments (nargs='?' so empty flags are allowed)
    p.add_argument(
        "-N", "--normal_bam",
        nargs='?', const='', default='',
        help="Input BAM (normal)"
    )
    p.add_argument(
        "--tumor_sample",
        nargs='?', const='', default='',
        help="Tumor sample name"
    )
    p.add_argument(
        "--normal_sample",
        nargs='?', const='', default='',
        help="Normal sample name"
    )
    p.add_argument(
        "--germline_resource",
        nargs='?', const='', default='',
        help="Germline resource VCF"
    )
    p.add_argument(
        "--panel_of_normals",
        nargs='?', const='', default='',
        help="Panel of normals VCF"
    )
    p.add_argument(
        "-L", "--intervals",
        nargs='?', const='', default='',
        help="Intervals/BED file"
    )
    p.add_argument(
        "--f1r2_tar_gz",
        nargs='?', const='', default='',
        help="F1R2 counts tar.gz output"
    )
    p.add_argument(
        "--bam_output",
        nargs='?', const='', default='',
        help="Realigned BAM output"
    )
    p.add_argument(
        "--native_pair_hmm_threads",
        nargs='?', type=int, const=0, default=0,
        help="Threads for native PairHMM"
    )
    return p.parse_args()



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
    if has(args.intervals):
        cmd += ["-L", args.intervals]
    if has(args.f1r2_tar_gz):
        cmd += ["--f1r2-tar-gz", args.f1r2_tar_gz]
    if has(args.bam_output):
        cmd += ["--bam-output", args.bam_output]
    if args.native_pair_hmm_threads and args.native_pair_hmm_threads > 0:
        cmd += ["--native-pair-hmm-threads", str(args.native_pair_hmm_threads)]
    cmd += ["-O", args.output_vcf]
    return cmd
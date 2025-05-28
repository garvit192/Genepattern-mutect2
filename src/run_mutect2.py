import subprocess
import sys

from preprocess import ensure_fai, ensure_dict, ensure_bam_index, ensure_feature_index, has
from cmd_builder import build_mutect2_cmd, parse_args


def main():
    args = parse_args()

    # always preprocess reference and tumor BAM
    ensure_fai(args.reference)
    ensure_dict(args.reference)
    ensure_bam_index(args.input_bam)

    # optional preprocessing
    if has(args.normal_bam):
        ensure_bam_index(args.normal_bam)
    if has(args.germline_resource):
        ensure_feature_index(args.germline_resource)
    if has(args.panel_of_normals):
        ensure_feature_index(args.panel_of_normals)

    # build and run Mutect2 command
    cmd = build_mutect2_cmd(args)
    print(f"[+] Running:\n    {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Mutect2 failed (exit {e.returncode})", file=sys.stderr)
        sys.exit(e.returncode)
        

if __name__ == "__main__":
    main()

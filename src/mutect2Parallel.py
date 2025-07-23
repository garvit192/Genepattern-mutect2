import subprocess
import sys
import copy
from pathlib import Path

from preprocess import ensure_fai, ensure_dict, ensure_bam_index, ensure_feature_index, has
from cmd_builder import build_mutect2_cmd, parse_args


def run_for_chromosome(base_args, chromosome):
    print(f"[+] Preparing Mutect2 run for {chromosome}")
    args = copy.deepcopy(base_args)
    args.intervals = chromosome
    args.output_vcf = f"{chromosome}_somatic.vcf"

    cmd = build_mutect2_cmd(args)
    print(f"[+] Running for {chromosome}: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print(f"[+] Completed: {chromosome}")
    except subprocess.CalledProcessError as e:
        print(f"[!] Failed on {chromosome} (exit {e.returncode})")


def main():
    
    args = parse_args()

    ensure_fai(args.reference)
    ensure_dict(args.reference)
    ensure_bam_index(args.input_bam)
    if has(args.normal_bam):
        ensure_bam_index(args.normal_bam)
    if has(args.germline_resource):
        ensure_feature_index(args.germline_resource)
    if has(args.panel_of_normals):
        ensure_feature_index(args.panel_of_normals)

    if args.parallel:
        from concurrent.futures import ProcessPoolExecutor
        from functools import partial
        from time import time

        CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]
        MAX_PARALLEL = 5

        start_time = time()
        run_partial = partial(run_for_chromosome, args)
        with ProcessPoolExecutor(max_workers=MAX_PARALLEL) as executor:
            executor.map(run_partial, CHROMOSOMES)
        print(f"[✓] Total time taken: {time() - start_time:.2f} seconds")
        return

    # Serial mode

    cmd = build_mutect2_cmd(args)
    print(f"[+] Running:\n    {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Mutect2 failed (exit {e.returncode})", file=sys.stderr)
        sys.exit(e.returncode)


if __name__ == "__main__":
    main()

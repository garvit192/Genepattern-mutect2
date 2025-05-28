import os
import subprocess
import sys

def has(arg):
    """Return True if arg is neither None nor an empty string."""
    return arg is not None and arg != ""


def ensure_fai(ref):
    if not os.path.isfile(ref):
        sys.exit(f"ERROR: reference FASTA not found: {ref}")
    fai = ref + ".fai"
    if not os.path.exists(fai):
        print(f"[+] Creating FASTA index: {fai}")
        subprocess.run(["samtools", "faidx", ref], check=True)


def ensure_dict(ref):
    dict_file = os.path.splitext(ref)[0] + ".dict"
    if not os.path.exists(dict_file):
        print(f"[+] Creating sequence dictionary: {dict_file}")
        subprocess.run([
            "gatk", "CreateSequenceDictionary",
            "-R", ref,
            "-O", dict_file
        ], check=True)


def ensure_bam_index(bam):
    if not os.path.isfile(bam):
        sys.exit(f"ERROR: BAM not found: {bam}")
    bai = bam + ".bai"
    if not os.path.exists(bai):
        print(f"[+] Creating BAM index: {bai}")
        subprocess.run(["samtools", "index", bam], check=True)


def ensure_feature_index(path):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: feature file not found: {path}")
    idx = path + (".tbi" if path.endswith(".gz") else ".idx")
    if not os.path.exists(idx):
        print(f"[+] Indexing feature file: {path}")
        subprocess.run(["gatk", "IndexFeatureFile", "-I", path], check=True)


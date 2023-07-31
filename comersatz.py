#!/usr/bin/env python3

import argparse
import math
import random
import subprocess
import sys
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seed", help="Set random seed for read selection")

    parser.add_argument(
        "-n", "--size", type=int, help="Size of final output metagenome"
    )

    parser.add_argument(
        "-i",
        "--illumina",
        required=True,
        nargs=3,
        action="append",
        metavar=("FWD", "REV", "PROB"),
        help="Paths to a pair of reads and the proportion of the final data they should constitute",
    )

    parser.add_argument(
        "-o", "--output", required=True, type=Path, help="Output directory"
    )

    args = parser.parse_args()

    args.illumina = [
        [Path(fwd), Path(rev), float(prob)] for (fwd, rev, prob) in args.illumina
    ]

    return args


def main():
    args = arguments()

    construct_illumina_metagenome(args.illumina, args.size, args.output, args.seed)


def calculate_required_read_count(total: int, proportion_of_total: float) -> int:
    return math.floor(total * proportion_of_total) or 1


def sample_reads(reads: Path, required_reads: int, paired_read_names: list[str], seed: int = 11) -> str:
    approximate_total_reads = estimate_read_count(reads)

    proportion_required = min(1.0, (required_reads / approximate_total_reads) * 1.5)

    grep_cmd = ("seqkit", "grep", "-f", "-", reads)

    grepped_reads = subprocess.run(grep_cmd, capture_output=True, text=True, input="\n".join(paired_read_names))

    sample_cmd = (
        "seqkit",
        "sample",
        "--rand-seed",
        str(seed),
        "-p",
        str(proportion_required)
    )

    sampled = subprocess.run(sample_cmd, capture_output=True, text=True, input=grepped_reads)

    head_cmd = ("seqkit", "head", "-n", str(required_reads))

    selected_reads = subprocess.run(
        head_cmd, capture_output=True, text=True, input=sampled.stdout
    )

    return selected_reads.stdout


def shuffle_reads(reads: str, seed: int = 11) -> str:
    shuffle_cmd = ("seqkit", "shuffle", "--rand-seed", str(seed))

    shuffled = subprocess.run(shuffle_cmd, capture_output=True, text=True, input=reads)

    return shuffled.stdout

def get_paired_read_names(fwd: Path, rev: Path) -> list[str]:
    
    fwd_names = get_read_names(fwd)
    rev_names = get_read_names(rev)

    return sorted(set(fwd_names).intersection(set(rev_names)))

def get_read_names(reads: Path) -> list[str]:

    cmd = ("seqkit", "--name", "--only-id", reads)

    return subprocess.run(cmd, text=True, capture_output=True).stdout.splitlines()

def estimate_read_count(reads: Path, assumption_length: int = 150):
    size_bytes = reads.stat().st_size

    approx_n_reads = (size_bytes / assumption_length) // 2

    return approx_n_reads


def construct_illumina_metagenome(illumina_triplets, total_output_reads, outdir, seed):
    outdir.mkdir(parents=True, exist_ok=True)

    out_fwd, out_rev = (
        outdir.joinpath(x).with_suffix(".fastq") for x in ("fwd", "rev")
    )

    out_fwd_reads, out_rev_reads = [], []

    for fwd, rev, prob in illumina_triplets:
        desired_reads = calculate_required_read_count(total_output_reads, prob)

        paired_read_names = get_paired_read_names(fwd, rev)

        out_fwd_reads.append(sample_reads(fwd, desired_reads, seed))
        out_rev_reads.append(sample_reads(rev, desired_reads, seed))

    fwd_reads = shuffle_reads("\n".join(out_fwd_reads), seed)
    rev_reads = shuffle_reads("\n".join(out_rev_reads), seed)

    with out_fwd.open("w") as f, out_rev.open("w") as r:
        f.write(fwd_reads)
        r.write(rev_reads)


def rename_illumina_reads(read_triplet):
    # TODO add this later; not necessary for mvp
    pass


if __name__ == "__main__":
    main()

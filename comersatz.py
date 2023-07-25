#!/usr/bin/env python3

import argparse
import random
import sys
import subprocess
import math
from pathlib import Path

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("--seed", help="Set random seed for read selection")

    parser.add_argument("-n", "--size", type=int, help="Size of final output metagenome")

    parser.add_argument("-i", "--illumina", required = True, nargs=3, action="append", metavar=("FWD", "REV", "PROB"), help="Paths to a pair of reads and the proportion of the final data they should constitute")

    parser.add_argument("-o", "--output", required = True, type=Path, help = "Output directory")

    args = parser.parse_args()

    return args


def main():

    args = arguments()

    construct_illumina_metagenome(args.illumina, args.size, args.output, args.seed)


def calculate_required_read_count(total: int, proportion_of_total: float) -> int:
    
    return math.floor(total * proportion_of_total) or 1


def sample_reads(reads: Path, required_reads: int, seed: int = 11) -> str:
    
    shuffle_cmd = ("seqkit", "sample", "--rand-seed", str(seed), reads)
    
    shuffled = subprocess.run(shuffle_cmd, capture_output=True, text=True)

    head_cmd = ("seqkit", "head", "-n", required_reads)

    selected_reads = subprocess.run(head_cmd, capture_output=True, text=True, stdin=shuffled)

    return selected_reads


def construct_illumina_metagenome(illumina_triplets, total_output_reads, outdir, seed):

    out_fwd, out_rev = (outdir.joinpath(x).with_suffix(".fastq") for x in ("fwd", "rev"))

    for fwd, rev, prob in illumina_triplets:
        
        desired_reads = calculate_required_read_count(total_output_reads, prob)
        
        fwd_sample = sample_reads(fwd, desired_reads, seed)
        rev_sample = sample_reads(rev, desired_reads, seed)

        with out_fwd.open("a") as f, out_rev.open("a") as r:
            f.write(fwd_sample)
            r.write(rev_sample)
        
def rename_illumina_reads(read_triplet):
    # TODO add this later; not necessary for mvp
    pass


if __name__ == "__main__":
    main()
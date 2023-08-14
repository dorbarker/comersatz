#!/usr/bin/env python3

import argparse
import collections
import gzip
import itertools
import math
import mimetypes
import random
from pathlib import Path

from Bio import SeqIO


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


def estimate_read_count(reads: Path, assumption_length: int = 150):
    size_bytes = reads.stat().st_size

    approx_n_reads = (size_bytes / assumption_length) // 2

    return approx_n_reads


def determine_read_count(reads: Path):
    _open = gzip.open if mimetypes.guess_type(reads)[1] == "gzip" else open

    with _open(reads, "rt") as f:
        records = SeqIO.parse(f, "fastq")
        last = collections.deque(enumerate(records), maxlen=1)[0]
        idx, _ = last
    return idx


def determine_selections(read_count: int, desired_reads: int, seed: int):
    indices = list(range(read_count))
    random.seed(seed)
    selected_indices = set(random.sample(indices, desired_reads))
    return selected_indices


def select_reads(reads: Path, selected_indices: list[int]) -> list[SeqIO.SeqRecord]:
    _open = gzip.open if mimetypes.guess_type(reads)[1] == "gzip" else open

    with _open(reads, "rt") as f:
        indexed_reads = enumerate(SeqIO.parse(f, "fastq"))
        selected = list(
            filter(lambda idx_rec: idx_rec[0] in selected_indices, indexed_reads)
        )

    return selected


def shuffle_reads(*reads, seed: int) -> list[SeqIO.SeqRecord]:
    merged = list(itertools.chain.from_iterable(reads))
    random.seed(seed)
    random.shuffle(merged)
    return [seqrec for (idx, seqrec) in merged]


def construct_illumina_metagenome(illumina_triplets, total_output_reads, outdir, seed):
    outdir.mkdir(parents=True, exist_ok=True)

    out_fwd, out_rev = (
        outdir.joinpath(x).with_suffix(".fastq") for x in ("fwd", "rev")
    )

    out_fwd_reads, out_rev_reads = [], []

    for fwd, rev, prob in illumina_triplets:
        desired_reads = calculate_required_read_count(total_output_reads, prob)
        input_read_count = determine_read_count(fwd)  # should be identical for rev
        selected_indices = determine_selections(input_read_count, desired_reads, seed)

        fwd_reads = select_reads(fwd, selected_indices)
        rev_reads = select_reads(rev, selected_indices)

        out_fwd_reads.append(fwd_reads)
        out_rev_reads.append(rev_reads)

    shuffled_fwd_reads = shuffle_reads(*out_fwd_reads, seed=seed)
    shuffled_rev_reads = shuffle_reads(*out_rev_reads, seed=seed)

    with out_fwd.open("w") as f, out_rev.open("w") as r:
        SeqIO.write(shuffled_fwd_reads, f, "fastq")
        SeqIO.write(shuffled_rev_reads, r, "fastq")


def rename_illumina_reads(read_triplet):
    # TODO add this later; not necessary for mvp
    pass


if __name__ == "__main__":
    main()

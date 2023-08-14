#!/usr/bin/env python3

import argparse
import collections
import csv
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

    parser.add_argument(
        "--no-rename",
        action="store_true",
        help="Do not rename the output FASTQ files with faked headers",
    )

    args = parser.parse_args()

    args.illumina = [
        [Path(fwd), Path(rev), float(prob)] for (fwd, rev, prob) in args.illumina
    ]

    return args


def main():
    args = arguments()

    construct_illumina_metagenome(
        args.illumina, args.size, args.output, not args.no_rename, args.seed
    )


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


def fake_illumina_name(
    record: SeqIO.SeqRecord, idx: int, read_number: int
) -> SeqIO.SeqRecord:
    assert read_number in (1, 2)

    fake_name = (
        f"SIM:001:112358:{idx}:{idx:05}:{idx:05}:{idx:05} {read_number}:N:0:GATTACA"
    )

    old_id = record.id

    record.id = record.name = fake_name
    record.description = ""

    return old_id, fake_name, record


def select_illumina_reads(
    reads: Path, desired_reads: int, read_count: int, seed: int
) -> list[SeqIO.SeqRecord]:
    indices = determine_selections(read_count, desired_reads, seed)

    selected_reads = select_reads(reads, indices)

    return selected_reads


def postprocess_illumina_reads(
    selected_reads: list[SeqIO.SeqRecord], orientation: str, rename: bool, seed: int
) -> list[SeqIO.SeqRecord]:
    if orientation == "fwd":
        orientation_num = 1
    elif orientation == "rev":
        orientation_num = 2
    else:
        raise ValueError(f"{orientation} is not a valid read orientation")

    shuffled = shuffle_reads(*selected_reads, seed=seed)
    if rename:
        oldnames, newnames, seqrecords = zip(
            *[
                fake_illumina_name(record, idx, orientation_num)
                for idx, record in enumerate(shuffled, 1)
            ]
        )
        lookup = list(zip(oldnames, newnames))
        return seqrecords, lookup
    else:
        lookup = []
        return shuffled, lookup


def construct_illumina_metagenome(
    illumina_triplets, total_output_reads, outdir, rename, seed
):
    outdir.mkdir(parents=True, exist_ok=True)

    out_fwd, out_rev = (
        outdir.joinpath(x).with_suffix(".fastq") for x in ("fwd", "rev")
    )

    out_fwd_reads, out_rev_reads = [], []

    for fwd, rev, prob in illumina_triplets:
        desired_reads = calculate_required_read_count(total_output_reads, prob)
        input_read_count = determine_read_count(fwd)  # should be identical for rev

        out_fwd_reads.append(
            select_illumina_reads(fwd, desired_reads, input_read_count, seed)
        )
        out_rev_reads.append(
            select_illumina_reads(rev, desired_reads, input_read_count, seed)
        )

    fwd_reads, lookup = postprocess_illumina_reads(out_fwd_reads, "fwd", rename, seed)
    rev_reads, lookup = postprocess_illumina_reads(out_rev_reads, "rev", rename, seed)

    with out_fwd.open("w") as f, out_rev.open("w") as r:
        SeqIO.write(fwd_reads, f, "fastq")
        SeqIO.write(rev_reads, r, "fastq")

    if lookup:
        with outdir.joinpath("fastq_lookup.tsv").open("w") as lookup_file:
            writer = csv.writer(lookup_file, delimiter="\t")
            for row in lookup:
                writer.writerow(row)


if __name__ == "__main__":
    main()

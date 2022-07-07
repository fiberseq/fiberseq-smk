#!/usr/bin/env python3
import pysam
import numpy as np
import argparse
import logging
import sys
import array


def get_mods_from_rec(rec, mods=[("A", 0, "a"), ("T", 1, "a")], mask=True):
    seq = np.frombuffer(bytes(rec.query_sequence, "utf-8"), dtype="S1")
    positions = []
    for mod in mods:
        if mod in rec.modified_bases:
            pos = np.array(rec.modified_bases[mod])[:, 0]
            positions.append(pos)
    if len(positions) < 1:
        return None, None, None
    methylated_positions = np.concatenate(positions)
    methylated_positions.sort(kind="mergesort")

    AT_mask = (seq == b"A") | (seq == b"T")
    AT_positions = np.argwhere(AT_mask).transpose()[0]

    binary = np.zeros(shape=len(seq), dtype=np.uint8)
    binary[methylated_positions] = 1

    if mask:
        # TODO check for other mods
        return binary[AT_mask], AT_positions, methylated_positions
    return binary, AT_positions, methylated_positions


def get_nucleosome_starts(rec):
    return np.array(rec.get_tag("ns"))


def get_nuceosome_lengths(rec):
    return np.array(rec.get_tag("nl"))


def get_nucleosomes(rec):
    return get_nucleosome_starts(rec), get_nuceosome_lengths(rec)


def parse():
    """Console script for fibertools."""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("bam", help="aligned bam file from actc")
    parser.add_argument(
        "-n",
        "--num-train",
        help="Number of fibers used to train HMM.",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "-o",
        "--out",
        help="Output bam or json file.",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    parser.add_argument("-t", "--threads", help="n threads to use", type=int, default=4)
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    return args


def extract(bam, args):
    for rec in bam.fetch(until_eof=True):
        ns, nl = get_nucleosomes(rec)


def main():
    args = parse()
    bam = pysam.AlignmentFile(args.bam, threads=args.threads, check_sq=False)
    extract(bam, args)


if __name__ == "__main__":
    main()

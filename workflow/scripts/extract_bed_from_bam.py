#!/usr/bin/env python3
import pysam
import numpy as np
import argparse
import logging
import sys
from numba import njit

# C+m
CPG_MODS = ["C", 0, "m"]
M6A_MODS = [("A", 0, "a"), ("T", 1, "a")]


def get_mod_pos_from_rec(rec, mods=M6A_MODS):
    positions = []
    for mod in mods:
        if mod in rec.modified_bases_forward:
            pos = np.array(rec.modified_bases_forward[mod])[:, 0]
            positions.append(pos)
    if len(positions) < 1:
        return None
    mod_positions = np.concatenate(positions)
    mod_positions.sort(kind="mergesort")
    # print(mod_positions)
    return mod_positions


def get_start_length_tags(rec, start_tag="ns", length_tag="nl"):
    if not rec.has_tag(start_tag) or not rec.has_tag(length_tag):
        return None, None
    starts = np.array(rec.get_tag(start_tag))
    lengths = np.array(rec.get_tag(length_tag))
    return starts, lengths


def get_nucleosomes(rec):
    return get_start_length_tags(rec, start_tag="ns", length_tag="nl")


def get_accessible(rec):
    return get_start_length_tags(rec, start_tag="as", length_tag="al")


# len = 7
# 0|1|2|3|4|5|6
# 6|5|4|3|2|1|0
# rev comp [1,5) = [2,6)
# new_st = len(7) - old_en(5)
# new_en = len(7) - old_st(1)
# @njit
def liftover_helper(aligned_pairs, sts, ens, query_length, is_reverse=False):
    read_pos = aligned_pairs[:, 0]
    ref_pos = aligned_pairs[:, 1]

    # flip positison if reversed
    if is_reverse:
        new_st = query_length - ens
        new_en = query_length - sts
        sts = new_st
        ens = new_en

    # TODO filter for sts and ends within alignment
    st_idxs = np.searchsorted(read_pos, sts, side="left")
    ed_idxs = np.searchsorted(read_pos, ens, side="left")
    ref_sts = ref_pos[st_idxs]
    ref_ens = ref_pos[ed_idxs]
    print(aligned_pairs)
    return ref_sts, ref_ens


def liftover(rec, sts, ens, aligned_pairs=None):
    if aligned_pairs is None:
        aligned_pairs = np.array(rec.get_aligned_pairs(matches_only=True))
    return liftover_helper(
        aligned_pairs,
        sts,
        ens,
        rec.query_length,
        is_reverse=rec.is_reverse,
    )


def write_bed12(rec, starts, output, lengths=None, aligned_pairs=None):
    if starts is None:
        return
    strand = "-" if rec.is_reverse else "+"
    passes = round(rec.get_tag("ec"))
    rgb = "0,0,0"

    # bed 6
    output.write(
        f"{rec.query_name}\t0\t{rec.query_length}\t{rec.query_name}\t{passes}\t{strand}\t"
    )
    if lengths is None:
        lengths = np.ones(starts.shape, dtype=int)
    bc = starts.shape[0]
    bs = ",".join(starts.astype(str))
    bl = ",".join(lengths.astype(str))
    # bed 12
    output.write(f"0\t{rec.query_length}\t{rgb}\t{bc}\t{bl}\t{bs}\n")


def extract(bam, args):
    for rec in bam.fetch(until_eof=True):
        aligned_pairs = None
        if args.reference:
            aligned_pairs = np.array(rec.get_aligned_pairs(matches_only=True))

        if args.nuc is not None:
            ns, nl = get_nucleosomes(rec)
            if ns is None:
                continue
            liftover(rec, ns, ns + nl, aligned_pairs=aligned_pairs)
            write_bed12(rec, ns, args.nuc, lengths=nl, aligned_pairs=aligned_pairs)
            break
        if args.acc is not None:
            acc_s, acc_l = get_accessible(rec)
            write_bed12(
                rec, acc_s, args.acc, lengths=acc_l, aligned_pairs=aligned_pairs
            )
        if args.m6a is not None:
            mod_pos = get_mod_pos_from_rec(rec, mods=M6A_MODS)
            write_bed12(rec, mod_pos, args.m6a, aligned_pairs=aligned_pairs)
        if args.cpg is not None:
            mod_pos = get_mod_pos_from_rec(rec, mods=CPG_MODS)
            write_bed12(rec, mod_pos, args.cpg, aligned_pairs=aligned_pairs)


def parse():
    """Console script for fibertools."""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("bam", help="(un)aligned bam file.")
    parser.add_argument(
        "-m",
        "--m6a",
        help="Output m6a bed12.",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "-c",
        "--cpg",
        help="Output m6a bed12.",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "-n",
        "--nuc",
        help="Output nucleosomes in bed12.",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "-a",
        "--acc",
        help="Output accessible stretches in bed12.",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Make output bed12 files use reference coordinates.",
        action="store_true",
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


def main():
    args = parse()
    bam = pysam.AlignmentFile(args.bam, threads=args.threads, check_sq=False)
    extract(bam, args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from tracemalloc import start
import pysam
import numpy as np
import argparse
import logging
from numba import njit
import tqdm

# C+m
CPG_MODS = [("C", 0, "m")]
M6A_MODS = [("A", 0, "a"), ("T", 1, "a")]
D_TYPE = np.int64


def get_mod_pos_from_rec(rec, mods=M6A_MODS):
    if rec.modified_bases_forward is None:
        return None
    positions = []
    for mod in mods:
        if mod in rec.modified_bases_forward:
            pos = np.array(rec.modified_bases_forward[mod], dtype=D_TYPE)[:, 0]
            positions.append(pos)
    if len(positions) < 1:
        return None
    mod_positions = np.concatenate(positions, dtype=D_TYPE)
    mod_positions.sort(kind="mergesort")
    return mod_positions


def get_start_length_tags(rec, start_tag="ns", length_tag="nl"):
    if not rec.has_tag(start_tag) or not rec.has_tag(length_tag):
        return None, None
    starts = np.array(rec.get_tag(start_tag), dtype=D_TYPE)
    lengths = np.array(rec.get_tag(length_tag), dtype=D_TYPE)
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
@njit
def liftover_helper(aligned_pairs, sts, ens, query_length, is_reverse=False):
    read_pos = aligned_pairs[:, 0]
    ref_pos = aligned_pairs[:, 1]

    # flip positison if reversed
    if is_reverse:
        new_st = query_length - ens
        new_en = query_length - sts
        sts = new_st[::-1]
        ens = new_en[::-1]

    # search of closest matching index
    st_idxs = np.searchsorted(read_pos, sts, side="left")
    en_idxs = np.searchsorted(read_pos, ens, side="left")
    ## remove things that are past the end of the read
    st_idxs[st_idxs >= ref_pos.shape[0]] = ref_pos.shape[0] - 1
    en_idxs[en_idxs >= ref_pos.shape[0]] = ref_pos.shape[0] - 1
    # get the ref positions
    ref_sts = ref_pos[st_idxs]
    ref_ens = ref_pos[en_idxs]

    # remove zero length liftovers (past start or end of alignment)
    keep_idx = ref_ens - ref_sts > 0
    # if (~keep_idx).any():
    #    logging.debug(f"removing {(~keep_idx).sum()} zero length liftovers.")
    return ref_sts[keep_idx], ref_ens[keep_idx]


def liftover(rec, sts, ens, aligned_pairs=None):
    if aligned_pairs is None:
        aligned_pairs = np.array(rec.get_aligned_pairs(matches_only=True), dtype=D_TYPE)
    return liftover_helper(
        aligned_pairs,
        sts,
        ens,
        rec.query_length,
        is_reverse=rec.is_reverse,
    )


@njit
def make_bed_blocks(starts, lengths, st, en):
    o_starts = ""
    o_lengths = ""
    bc = 0
    if starts[0] != 0:
        o_starts += "0,"
        o_lengths += "1,"
        bc += 1
    for b_st, b_ln in zip(starts, lengths):
        o_starts += str(b_st) + ","
        o_lengths += str(b_ln) + ","
        bc += 1
    if starts[-1] + lengths[-1] != en:
        o_starts += str(en - st - 1) + ","
        o_lengths += "1,"
        bc += 1
    return bc, o_starts[:-1], o_lengths[:-1]


def write_bed12(rec, starts, output, lengths=None, aligned_pairs=None):
    if starts is None:
        return
    strand = "-" if rec.is_reverse else "+"
    passes = round(rec.get_tag("ec"))
    rgb = "0,0,0"

    if aligned_pairs is None:
        ct, st, en = rec.query_name, np.int64(0), np.int64(rec.query_length)
    else:
        ct, st, en = (
            rec.reference_name,
            np.int64(rec.reference_start),
            np.int64(rec.reference_end),
        )
    # add lengths if none are provided
    if lengths is None:
        lengths = np.ones(starts.shape, dtype=D_TYPE)

    if aligned_pairs is not None:
        l_sts, l_ens = liftover(
            rec, starts, starts + lengths, aligned_pairs=aligned_pairs
        )
        l_len = l_ens - l_sts
        l_sts -= st
        starts = l_sts
        lengths = l_len

    assert (
        starts + lengths > en
    ).sum() == 0, f"Blocks exceed end position of bed12\n{starts}\n{lengths}\n{en}"

    bc, bs, bl = make_bed_blocks(starts, lengths, st, en)
    if bs[0] != "0":
        logging.warning("First block start is not 0")
        return

    # write bed 6
    output.write(f"{ct}\t{st}\t{en}\t{rec.query_name}\t{passes}\t{strand}\t")
    # think start and end
    output.write(f"{st}\t{en}\t{rgb}\t")
    # bed 12
    output.write(f"{bc}\t{bl}\t{bs}\n")


def extract(bam, args):
    for rec in tqdm.tqdm(bam.fetch(until_eof=True)):
        aligned_pairs = None
        if args.reference and rec.is_unmapped:
            continue
        if args.reference:
            aligned_pairs = np.array(
                rec.get_aligned_pairs(matches_only=True), dtype=D_TYPE
            )

        if args.nuc is not None:
            ns, nl = get_nucleosomes(rec)
            write_bed12(rec, ns, args.nuc, lengths=nl, aligned_pairs=aligned_pairs)
            # break
        if args.msp is not None:
            msp_s, msp_l = get_accessible(rec)
            write_bed12(
                rec, msp_s, args.msp, lengths=msp_l, aligned_pairs=aligned_pairs
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
        help="Output CpG bed12.",
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
        "--msp",
        help="Output accessible stretches (MSPs) in bed12.",
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

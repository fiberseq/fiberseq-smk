#!/usr/bin/env python3
from operator import mod
import pysam
import sys
import numpy as np
from tqdm import tqdm
import pandas as pd
import numpy as np
import argparse
import logging
import pomegranate as pom
import tqdm


def get_mods_from_rec(rec, mods=[("A", 0, "a"), ("T", 1, "a")], mask=True, binary=True):
    seq = np.frombuffer(bytes(rec.query_sequence, "utf-8"), dtype="S1")
    positions = []
    for mod in mods:
        if mod in rec.modified_bases:
            pos = np.array(rec.modified_bases[mod])[:, 0]
            positions.append(pos)
    methylated_positions = np.concatenate(positions)
    methylated_positions.sort(kind="mergesort")
    if binary:
        binary = np.zeros(shape=len(seq), dtype=bool)
        binary[methylated_positions] = 1
        if mask:
            # TODO check for other mods
            return binary[(seq == b"A") | (seq == b"T")]
        return binary
    return methylated_positions


def train_hmm(data, output):
    model = pom.HiddenMarkovModel().from_samples(
        pom.DiscreteDistribution,
        n_components=2,
        X=data,
        verbose=True,
        n_jobs=1,
        max_iterations=10,  # 250,
        n_init=10,
    )
    model.bake()
    return model


def parse():
    """Console script for fibertools."""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("bam", help="aligned bam file from actc")
    parser.add_argument(
        "-m",
        "--model",
        help="pretrained hmm model (json).",
        default=None,
    )
    parser.add_argument(
        "-n",
        "--num-train",
        help="Number of fibers used to train HMM.",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "out",
        help="Output bam or json file.",
        type=argparse.FileType("w"),
    )
    parser.add_argument("-t", "--threads", help="n threads to use", type=int, default=1)
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
    if args.model is None:
        training_set = []
        for rec in bam.fetch(until_eof=True):
            mods = get_mods_from_rec(rec, binary=True, mask=True)
            training_set.append(mods)
        model = train_hmm(training_set, args.out)
        json_model = model.to_json()
        with args.out as handle:
            handle.write(json_model)
    else:
        out = pysam.AlignmentFile(args.out, "wb", template=bam)

    return 0


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import pysam
import numpy as np
import argparse
import logging
import pomegranate as pom
from hmmutils import *


def train_hmm(data, n_jobs=4):
    model = pom.HiddenMarkovModel().from_samples(
        pom.DiscreteDistribution,
        n_components=2,
        X=data,
        verbose=True,
        max_iterations=250,
        n_init=10,
        n_jobs=n_jobs,
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
        "-c", "--cutoff", type=int, default=65, help="hmm nucleosome size cutoff"
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
        for idx, rec in enumerate(bam.fetch(until_eof=True)):
            mods, _AT_pos, _m6a_pos = get_mods_from_rec(rec, mask=True)
            training_set.append(mods)
            if idx >= args.num_train:
                break
        model = train_hmm(training_set, n_jobs=args.threads)
        json_model = model.to_json()
        with args.out as handle:
            handle.write(json_model)
    else:
        out = pysam.AlignmentFile(args.out, "wb", template=bam)
        hmm = pom.HiddenMarkovModel().from_json(args.model)
        _actuated_label, nucleated_label = assign_states(hmm)
        apply_hmm(bam, hmm, nucleated_label, args.cutoff, out)

    return 0


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from distutils.log import debug
from statistics import mode
import pysam
import sys
import numpy as np
from tqdm import tqdm
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import argparse
import logging
import tqdm
import pickle


def coordinateConversion_MMTag(sequence, base, modification_coords):
    """Sequence is array of bases. Base is a base
    for conversion to the coordinate system
    used in the MM tag."""

    mask = sequence == bytes(base, encoding="utf-8")
    # find all masks = boolean Array
    coords = modification_coords[
        sequence[modification_coords] == bytes(base, encoding="utf-8")
    ]
    # when working with double stranded data we can only use modifications
    # on the specifc base that the modification would fall on on that strand
    # i.e. As  on + , Ts on -, we only want the mods that = that base of interest
    masked_coords = np.cumsum(mask)[coords] - 1
    distance = np.diff(masked_coords) - 1
    if coords.shape[0] < 1:
        return "", 0
    return (
        ",".join([masked_coords[0].astype(str)] + list(distance.astype(str))),
        coords.shape[0],
    )


def train_gmm(ipdRatios, means_init=[0.5, 2.2]):
    gmm = GaussianMixture(
        n_components=2,
        n_init=3,
        max_iter=500,
        covariance_type="full",
        tol=1e-5,
        means_init=means_init,
    )
    return gmm.fit(ipdRatios)


def apply_model(model, tpl, ipdRatios, min_prediction_value):
    labels = model.predict_proba(ipdRatios)
    predictions = labels[:, np.argmax(model.means_[:, 0])]
    return tpl[predictions >= min_prediction_value]


def write_model(model, file):
    with open(file, "wb") as f:
        pickle.dump(model, f)
    logging.debug("Done writing model")


def read_model(filename):
    return pickle.load(open(filename, "rb"))


def apply_gmm(
    csv,
    bam,
    out,
    min_prediction_value=0.99999999,
    min_number_of_calls=25,
    pre_trained_model=None,
):
    for rec in tqdm.tqdm(bam.fetch(until_eof=True), total=csv.index.unique().shape[0]):
        if not rec.query_name in csv.index:
            logging.debug(f"Missing {rec.query_name}")
            out.write(rec)
            continue
        molecule_df = csv.loc[rec.query_name]

        # if there are not enough m6a calls, skip
        if molecule_df.shape[0] < min_number_of_calls:
            out.write(rec)
            continue

        # convert to zero based coordinates on the read
        tpl = molecule_df.tpl.to_numpy() - 1
        # ipd ratios for the above positions
        ipdRatios = molecule_df.ipdRatio.to_numpy().reshape(-1, 1)

        # either do per molecule training or use a pre-trained model
        if pre_trained_model is None:
            model = train_gmm(ipdRatios)
        else:
            logging.debug("Using pre-trained model")
            model = pre_trained_model

        # coordinates of m6a calls predicted by the gmm
        mol_m6a = apply_model(model, tpl, ipdRatios, min_prediction_value)

        # check the mod count is correct
        sequence = np.frombuffer(bytes(rec.query_sequence, "utf-8"), dtype="S1")
        mod_count = ((sequence[mol_m6a] == b"A") | (sequence[mol_m6a] == b"T")).sum()
        assert mol_m6a.shape[0] == mod_count, f"{mol_m6a.shape[0]} != {mod_count}"
        if mod_count < 1:
            out.write(rec)
            continue

        # add modifications to bam
        A_mods, A_mod_count = coordinateConversion_MMTag(sequence, "A", mol_m6a)
        T_mods, T_mod_count = coordinateConversion_MMTag(sequence, "T", mol_m6a)
        logging.debug(
            f"{rec.query_name} has {A_mod_count:,} A mods and {T_mod_count:,} T mods"
        )
        mods = ""
        if A_mod_count > 0:
            mods += "A+a," + A_mods + ";"
        if T_mod_count > 0:
            mods += "T-a," + T_mods + ";"
        new_probabilities = [255] * mod_count

        assert (
            A_mod_count + T_mod_count == mod_count
        ), f"{A_mod_count + T_mod_count} != {mod_count}"

        assert (
            len(new_probabilities) == mod_count
        ), f"{len(new_probabilities)} != {mod_count}"

        # check if tag exists
        if rec.has_tag("MM"):
            original_mods = rec.get_tag("MM")
            new_mods = str(original_mods) + mods

            original_probs = list(rec.get_tag("ML"))
            new_probs = (
                np.array(original_probs + new_probabilities).astype(int).tolist()
            )

            new_tags = [("MM", new_mods, "Z"), ("ML", new_probs, "C")]
            for tag_set in rec.get_tags():
                if tag_set[0] != "MM" and tag_set[0] != "ML":
                    new_tags.append(tag_set)
            rec.set_tags(new_tags)
        else:
            rec.tags += [("MM:Z:", mods)]
            rec.tags += [("ML:C:", new_probabilities)]
        out.write(rec)


def read_csv(file):
    csv = pd.read_csv(
        file,
        usecols=[0, 1, 2, 3, 8, 9],
        dtype={
            "refName": str,
            "tpl": int,
            "strand": bool,
            "base": str,
            "ipdRatio": float,
            "coverage": int,
        },
    )
    csv = csv[csv.base == "A"]
    csv.set_index("refName", inplace=True)
    logging.debug("Done reading csv")
    return csv


def parse():
    """Console script for fibertools."""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("csv", help="csv file with tpl, ipdRatio, and coverage")
    parser.add_argument("bam", help="aligned bam file from actc")
    parser.add_argument(
        "-p",
        "--min-prediction-value",
        help="Minimum prediction value from the gmm to set m6a.",
        type=float,
        default=0.99999999,
    )
    parser.add_argument(
        "-m",
        "--min-number-of-calls",
        help="Min number of calls per fiber to use.",
        type=int,
        default=25,
    )
    parser.add_argument(
        "--model",
        help="Pre trained GMM to load.",
        default=None,
    )
    parser.add_argument(
        "--train",
        help="train a GMM instead of writing a bam",
        action="store_true",
    )
    parser.add_argument("-o", "--out", help="Output bam file.", default=sys.stdout)
    parser.add_argument("-t", "--threads", help="n threads to use", type=int, default=1)
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    return args


def train_model_on_csv(csv, args):
    ipdRatios = csv.ipdRatio.to_numpy().reshape(-1, 1)
    logging.debug(f"Starting training on {ipdRatios.shape[0]:,} m6a calls")
    model = train_gmm(ipdRatios)
    write_model(model, args.out)


def main():
    args = parse()
    csv = read_csv(args.csv)

    if args.model is not None:
        logging.debug(f"Loading pre-trained model from {args.model}")
        model = read_model(args.model)
    else:
        model = None

    if args.train:
        train_model_on_csv(csv, args)
    else:
        bam = pysam.AlignmentFile(args.bam, threads=args.threads, check_sq=False)
        out = pysam.AlignmentFile(args.out, "wb", template=bam)
        apply_gmm(
            csv,
            bam,
            out,
            min_prediction_value=args.min_prediction_value,
            min_number_of_calls=args.min_number_of_calls,
            pre_trained_model=model,
        )
    return 0


if __name__ == "__main__":
    main()

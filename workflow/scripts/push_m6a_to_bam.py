#!/usr/bin/env python3
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

    MM_coords = ",".join(list((np.diff(np.cumsum(mask)[coords]) - 1).astype(str)))

    return MM_coords


def apply_gmm(csv, bam, out, min_prediction_value=0.99999999, min_number_of_calls=25):
    for rec in tqdm.tqdm(bam.fetch(until_eof=True), total=csv.index.unique().shape[0]):
        if not rec.query_name in csv.index:
            logging.debug(f"Missing {rec.query_name}")
            out.write(rec)
            continue
        molecule_df = csv.loc[rec.query_name]
        # ec = molecule_df["coverage"].mean()  # effective coverage
        # data = molecule_df[["tpl", "ipdRatio"]].to_numpy().transpose()
        # tpl = data[0].astype(int)
        # ipdRatios = data[1].reshape(-1, 1)

        # convert to zero based
        tpl = molecule_df.tpl.to_numpy() - 1
        ipdRatios = molecule_df.ipdRatio.to_numpy().reshape(-1, 1)
        # logging.debug(f"{rec.query_name} {tpl.shape} {ipdRatios.shape}")
        if int(ipdRatios.shape[0]) < min_number_of_calls:
            continue

        gmm = GaussianMixture(
            n_components=2, n_init=3, max_iter=500, covariance_type="full", tol=1e-5
        )
        model = gmm.fit(ipdRatios)

        labels = model.predict_proba(ipdRatios)

        predictions = labels[:, np.argmax(model.means_[:, 0])]
        mol_m6a = tpl[predictions >= min_prediction_value]

        # check the mod count is correct
        sequence = np.frombuffer(bytes(rec.query_sequence, "utf-8"), dtype="S1")
        mod_count = ((sequence[mol_m6a] == b"A") | (sequence[mol_m6a] == b"T")).sum()
        assert mol_m6a.shape[0] == mod_count, f"{mol_m6a.shape[0]} != {mod_count}"

        # add modifications to bam
        A_mods = coordinateConversion_MMTag(sequence, "A", mol_m6a)
        T_mods = coordinateConversion_MMTag(sequence, "T", mol_m6a)
        mods = "A+a," + A_mods + ";" + "T-a," + T_mods + ";"
        new_probabilities = [255] * mod_count
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
            rec.tags += [("ML:Z:", new_probabilities)]
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


def main():
    args = parse()
    csv = read_csv(args.csv)
    bam = pysam.AlignmentFile(args.bam, threads=args.threads, check_sq=False)
    out = pysam.AlignmentFile(args.out, "wb", template=bam)
    apply_gmm(
        csv,
        bam,
        out,
        min_prediction_value=args.min_prediction_value,
        min_number_of_calls=args.min_number_of_calls,
    )
    return 0


if __name__ == "__main__":
    main()
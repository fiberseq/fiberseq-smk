# build pipeline to generate per fiber GMM model
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import normalize
import sys
import time
import tqdm
import pysam
import multiprocessing as mp

# line format
# m64018_201129_132425/18/ccs",2,0,G,6,0.814,0.248,0.639,1.273,11

# first arg = csv
# second arg = aligned bam
# third arg = file_name
# fourth arg = number of cpus to use

# save dictionary of tpl_arrays per molecule name to then
# use in converting to genomic coordinates
# but first write out bed file to file (unaligned fiber coordinates)


def trainGMMs(csv):
    intervals = []
    molecule_dict = {}
    # 	print('here')
    for molecule in tqdm.tqdm(csv.groupby("refName")):
        # 		print(molecule)
        molecule_name = molecule[0]
        molecule_df = molecule[1]

        ec = np.mean(molecule_df["coverage"])  # effective coverage

        data = molecule_df[["tpl", "ipdRatio"]].to_numpy().transpose()

        tpl = data[0].astype(int)
        ipdRatios = data[1].reshape(-1, 1)
        if int(ipdRatios.shape[0]) < 25:
            continue

        gmm = GaussianMixture(
            n_components=2, n_init=3, max_iter=500, covariance_type="full", tol=1e-5
        )
        model = gmm.fit(ipdRatios)

        labels = model.predict_proba(ipdRatios)

        predictions = labels[:, np.argmax(model.means_[:, 0])]
        fiber_positions_where_methylated = tpl[predictions >= 0.99999999]

        # full_name = str(molecule_df.refName.iloc[0])
        molecule_dict[
            int(molecule_name.split("/")[1])
        ] = fiber_positions_where_methylated

        chrom = molecule_name

        start, stop = min(tpl), max(tpl)

        blockCount = str(len(fiber_positions_where_methylated) + 2)

        blockStarts = ",".join(
            ["0"]
            + list((fiber_positions_where_methylated - 1).astype(str))
            + [str((stop - start) - 1)]
        )

        blockSizes = ",".join(list(np.ones(int(blockCount), dtype=int).astype(str)))

        zmwid = str(molecule_name)

        bed_interval = [
            chrom,  # chromosome
            str(start),  # fiber start in genomic coords
            str(stop),  # fiber stop in genomic coords
            zmwid,  # ZMW ID
            str(ec),  # effective coverage
            ".",
            str(start),  # fiber start in genomic coords
            str(stop),  # fiber stop in genomic coords]
            "128,0,128",
            blockCount,
            blockSizes,
            blockStarts,
        ]

        intervals.append("\t".join(bed_interval))

    return molecule_dict, intervals


def main():

    prefix = sys.argv[3]

    ### first read in just zmwids
    ### then read in csv for each block
    ### avoiding the unused zmwids
    ### and use multiprocessing
    ### merge dictionaries
    ### and merge lists

    csv = pd.read_csv(
        sys.argv[1],
        usecols=[0, 1, 3, 8, 9],
        names=["refName", "tpl", "base", "ipdRatio", "coverage"],
        dtype={
            "refName": str,
            "tpl": int,
            "base": str,
            "ipdRatio": float,
            "coverage": int,
        },
    )
    csv = csv[csv.base == "A"]

    molecule_to_zmwid, molecule_file = trainGMMs(csv)

    molecule_out_name = prefix + ".unaligned.bed"
    with open(molecule_out_name, "w") as handle:
        handle.write("\n".join(molecule_file))
        handle.write("\n")

    # generate aligned coordinates

    aligned_file = []
    # iterate through bam and generate genome specific bed file
    bam = pysam.AlignmentFile(sys.argv[2], "rb")

    for interval in tqdm.tqdm(bam):

        zmwid = int(interval.get_tag("zm"))

        if zmwid in molecule_to_zmwid:

            pos = molecule_to_zmwid[zmwid]

            aligned_pairs = (
                np.array(interval.get_aligned_pairs(matches_only=True))
                .astype(int)
                .transpose()
            )

            start, stop = int(interval.reference_start), int(interval.reference_end)

            if interval.is_reverse:

                # if the read is reversed. The coordinates from the csv
                # are in reverse (from the end of the read)
                # we need to change the positions we get fro the csv
                # file before we look at aligned pairs

                pos = np.abs(pos - int(interval.infer_read_length()))
                ref_coords = aligned_pairs[1][np.isin(aligned_pairs[0], pos)] - start

            else:

                ref_coords = (
                    aligned_pairs[1][np.isin(aligned_pairs[0], pos)] - start
                ) - 1

            chrom = str(interval.reference_name)

            smrt_cell = str(interval.query_name.split("/")[0])

            blockCount = str(len(ref_coords) + 2)

            blockStarts = ",".join(
                ["0"] + list((ref_coords).astype(str)) + [str((stop - start) - 1)]
            )

            blockSizes = ",".join(list(np.ones(int(blockCount), dtype=int).astype(str)))

            bed_interval = [
                chrom,  # chromosome
                str(start),  # fiber start in genomic coords
                str(stop),  # fiber stop in genomic coords
                str(zmwid) + "/" + smrt_cell,  # fiber name (zmwid + "/" + smrt cell )
                str(interval.get_tag("ec")),  # effective coverage
                ".",
                str(start),  # fiber start in genomic coords
                str(stop),  # fiber stop in genomic coords]
                "128,0,128",
                blockCount,
                blockSizes,
                blockStarts,
            ]

            aligned_file.append("\t".join(bed_interval))

    aligned_out_name = prefix + ".aligned.bed"

    with open(aligned_out_name, "w") as handle:
        handle.write("\n".join(aligned_file))
        handle.write("\n")

        # print('\t'.join(bed_interval))
        # print(ref_coords)


if __name__ == "__main__":
    main()

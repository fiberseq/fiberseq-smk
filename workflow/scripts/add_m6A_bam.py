import pysam
from pybedtools import BedTool
import sys
import numpy as np
from tqdm import tqdm
import pyfaidx as pf


bam = pysam.AlignmentFile(sys.argv[1], "rb")  # aligned bam
bed = BedTool(sys.argv[2])  # aligned bed

output_bam_name = sys.argv[1][:-3] + "m6A.tagged.bam"

header = bam.header.to_dict()

output_bam = pysam.AlignmentFile(output_bam_name, "wb", header=header)

bed_dict = {}

# bed to dictionary conversion for random access to methylation
# marks encoded in the bed file

for interval in bed:

    name = str(interval.name)

    methylations = np.fromstring(interval[-1], sep=",", dtype=np.float64)

    bed_dict[name] = (methylations, int(interval.start))


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


with output_bam as out_f:  # open output bam file

    with bam as in_f:  # open input bam file, iterate through,
        # and write modification to output bam
        for read in tqdm(in_f):

            name_list = read.query_name.split("/")
            new_name = str(name_list[1] + "/" + name_list[0])

            # new_name = str(read.query_name)

            methylations = bed_dict.get(new_name)

            if type(methylations) != type(None):

                ap = np.vstack(read.get_aligned_pairs(matches_only=True)).T

                # adjust bed to genomic coordinates
                # search for adjusted coordinates
                # in ap -- take the molecular coordinates that aligned

                genomic_m6a = methylations[0].astype(int) + methylations[1]

                mol_m6a = ap[
                    0, np.isin(ap[1], genomic_m6a)
                ]  # this is coordinates of query sequence
                # we need coordinates of query alignemnt sequence

                sequence = np.frombuffer(
                    bytes(read.query_sequence, "utf-8"), dtype="S1"
                )

                ##########
                # np.sum(np.isin(sequence[mol_m6a[1:-1]],['G','C']))
                # -- this will retrieve count of methylation calls with poor alignment
                # where reference does not match up A/T
                #########

                # # need to encode per strand
                # # A+a encodes forward strand ( check As , offset of As )
                # # T-a encodes reverse strand ( check Ts , offset of Ts )

                # need to grab all A and all T positions
                # and then frame methylation positions in
                # the appropriate base-space

                A_mods = coordinateConversion_MMTag(sequence, "A", mol_m6a[1:-1])
                T_mods = coordinateConversion_MMTag(sequence, "T", mol_m6a[1:-1])

                # need to check if MM tag exists
                # if "MM" tag exists -- we are going to
                # append to the existing tag
                # otherwise, we initialize the tag

                mod_count = A_mods.count(",") + 1 + T_mods.count(",") + 1

                mods = "A+a," + A_mods + ";" + "T-a," + T_mods + ";"
                new_probabilities = [255] * mod_count
                # check if tag exists

                if read.has_tag("MM"):
                    # if has tag
                    # we need to modify the flag

                    original_mods = read.get_tag("MM")
                    new_mods = str(original_mods) + mods

                    original_probs = list(read.get_tag("ML"))
                    # 					print(original_probs)
                    new_probs = (
                        np.array(original_probs + new_probabilities)
                        .astype(int)
                        .tolist()
                    )

                    # read.set_tag('MM', new_mods, value_type = 'Z', replace = True)
                    # read.set_tag('ML', new_probs, value_type = 'B', replace = True)

                    new_tags = [("MM", new_mods, "Z"), ("ML", new_probs, "C")]

                    for tag_set in read.get_tags():
                        if tag_set[0] != "MM" and tag_set[0] != "ML":

                            new_tags.append(tag_set)

                    read.set_tags(new_tags)

                else:
                    # if does not have tag just add it
                    read.tags += [("MM:Z:", mods)]

                out_f.write(read)

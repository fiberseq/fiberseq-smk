#!/usr/bin/env python3
import logging
import argparse

# https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n))


def parse():
    """Split zmw list into many files"""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("zmw", help="aligned bam file from actc")
    parser.add_argument(
        "-o",
        "--out",
        help="Output files",
        nargs="+",
    )
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
    zmws = [line for line in open(args.zmw, "r")]
    n_chunks = len(args.out)
    for zmw_batch, out_file in zip(split(zmws, n_chunks), args.out):
        with open(out_file, "w") as f:
            for zmw in zmw_batch:
                f.write(zmw)


if __name__ == "__main__":
    main()

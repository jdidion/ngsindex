import argparse
import ngsindex
import json
from xphyle import STDOUT, open_


def summarize():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--output", default=STDOUT,
        help="Output file (defaults to stdout)."
    )
    parser.add_argument("index", help="Index file.")
    args = parser.parse_args()

    index = ngsindex.parse_index(args.index)
    summary = index.summarize()
    with open_(args.output, "w") as out:
        json.dump(summary, out, indent=4)

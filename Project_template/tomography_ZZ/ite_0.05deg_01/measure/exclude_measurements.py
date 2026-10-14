#!/usr/bin/env python3
"""Exclude phase delay measurements listed in a droplist.

A measurement is identified by its source station, receiver station and
frequency tag. These are parsed from each delay-file line after replacing
'/' and '_' with spaces (fields 4, 5 and 12, 1-indexed). Each droplist line
is split the same way, and a measurement is dropped if some droplist line
contains all three values as exact fields (no substring matching).
"""
import argparse
import sys


def fields(line):
    """Split a line into fields, treating '/' and '_' as whitespace."""
    return line.replace("/", " ").replace("_", " ").split()


def parse_args():
    p = argparse.ArgumentParser(
        description="Exclude phase delay measurements (station pair + frequency tag) "
                    "that appear in a droplist."
    )
    p.add_argument("-i", "--input", required=True, help="input delay (measurement) file")
    p.add_argument("-d", "--droplist", required=True, help="droplist file")
    p.add_argument("-o", "--output", required=True, help="output delay file (droplist delays removed)")
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.droplist) as f:
        drops = [set(fields(ln)) for ln in f]

    kept = dropped = 0
    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            f = fields(line)
            if len(f) < 12:
                # Can't identify the measurement; keep it.
                fout.write(line + "\n")
                kept += 1
                continue

            key = {f[3], f[4], f[11]}  # source, receiver, frequency tag
            if any(key <= d for d in drops):
                dropped += 1
            else:
                fout.write(line + "\n")
                kept += 1

    print(f"kept {kept}, excluded {dropped} -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

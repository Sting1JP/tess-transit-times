#!/usr/bin/env python3
"""Compare two times.csv files written by timing.py, for example a rerun with a known answer.

    python3 tools/compare_times.py known_answers/wasp22/times.csv output/timing/wasp22/DATE_TIME/times.csv

(DATE_TIME is the run folder that timing.py prints at the end.) For the transits timed in both files
it prints the rms and the largest difference of the mid-times, and the largest difference of the
errors, in seconds, and it lists any transit timed in only one file. The exit status is 0 when both
files time the same transits, 1 when they do not or neither file has a timed transit, and 2 when a
file cannot be read as a times.csv file.
"""
import csv
import math
import sys

COLUMNS = ("epoch_tool", "t_mid_bjd_tdb", "err_d")


def load(path):
    times = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or any(c not in reader.fieldnames for c in COLUMNS):
            raise ValueError(f"{path}: no {', '.join(COLUMNS)} columns")
        for row in reader:
            if row["t_mid_bjd_tdb"]:
                times[int(row["epoch_tool"])] = (float(row["t_mid_bjd_tdb"]), float(row["err_d"]))
    return times


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    try:
        a, b = load(argv[0]), load(argv[1])
    except (OSError, KeyError, ValueError, TypeError, csv.Error) as exc:
        print(f"cannot read the files as timing.py times.csv files: {exc!r}")
        return 2
    if not a and not b:
        print("neither file has a timed transit")
        return 1
    common = sorted(set(a) & set(b))
    if common:
        dt = [(b[n][0] - a[n][0]) * 86400.0 for n in common]
        de = [(b[n][1] - a[n][1]) * 86400.0 for n in common]
        bad = [n for n, x, y in zip(common, dt, de) if not (math.isfinite(x) and math.isfinite(y))]
        dt = [x for x in dt if math.isfinite(x)]
        de = [y for y in de if math.isfinite(y)]
        print(f"{len(common)} transits timed in both files")
        if dt:
            rms = math.sqrt(sum(x * x for x in dt) / len(dt))
            print(f"mid-time, second minus first: rms {rms:.2f} s, largest {max(dt, key=abs):+.2f} s")
        if de:
            print(f"error, second minus first: largest {max(de, key=abs):+.2f} s")
        if bad:
            print("not a finite number in one of the files, epochs:", " ".join(map(str, bad)))
    only_a, only_b = sorted(set(a) - set(b)), sorted(set(b) - set(a))
    if only_a:
        print("timed only in the first file, epochs:", " ".join(map(str, only_a)))
    if only_b:
        print("timed only in the second file, epochs:", " ".join(map(str, only_b)))
    same = not only_a and not only_b
    print("the same transits are timed in both files" if same else "the files time different transits")
    return 0 if same else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

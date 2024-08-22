from datetime import datetime
import sys

sys.path.insert(0, '/home/sebschmi/git/short-range-template-switching/tsm_methods/tools/fpa')
import fpa_ext2

import argparse
parser = argparse.ArgumentParser(description="Benchmark fpa")
parser.add_argument("--length", type=int, default=200, help="Sequence length")

args = parser.parse_args()
length = args.length

lines = []
with open('sequences.txt') as file:
    lines = [line.strip()[:length] for line in file.readlines()]

for line in lines:
    assert len(line) == length, "sequence not long enough"

start_time = datetime.now()

f = fpa_ext2.FPA2()
f.set_int("min_length",8)
f.set_int("scan_flank",100)
f.set_int("scan_flank",60)

for [line1, line2] in zip(lines[1:], lines[:-1]):
    hits1 = f.scan_two(line1, line2, False)
    hits2 = f.scan_two(line2, line1, False)

end_time = datetime.now()

print(f"Computing {len(lines) - 1} alignments of length {length} took {(end_time - start_time).total_seconds()} seconds")

#!/bin/bash
#SBATCH --job-name=fix_ncbi
#SBATCH --output=fix_ncbi_%j.out
#SBATCH --error=fix_ncbi_%j.err
#SBATCH --cpus-per-task=1

set -euo pipefail

# Usage:
# sbatch fix_ncbi_fasta.slurm input.fasta output.fasta

INPUT="$1"
OUTPUT="$2"

python3 - "$INPUT" "$OUTPUT" <<'PY'
import sys

input_fa  = sys.argv[1]
output_fa = sys.argv[2]

# Entire contaminant contigs to remove
exclude = {
    "ctg_388_1"
}

# NCBI coordinates: 1-based, inclusive
# adaptor regions to remove by splitting the contig
split_regions = {
    "ctg_48_1": (1390, 1473),
    "ctg_7_1":  (17702306, 17702389)
}

MIN_LEN = 200


def read_fasta(path):
    header = None
    seq = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)

                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())

        if header is not None:
            yield header, "".join(seq)


def trim_terminal_ns(seq):
    original_len = len(seq)

    seq = seq.lstrip("Nn")
    left = original_len - len(seq)

    before_right = len(seq)
    seq = seq.rstrip("Nn")
    right = before_right - len(seq)

    return seq, left, right


def write_record(out, header, seq):
    out.write(f">{header}\n")
    for i in range(0, len(seq), 80):
        out.write(seq[i:i+80] + "\n")


input_records = 0
output_records = 0

input_bp = 0
output_bp = 0

excluded_bp = 0
adapter_bp = 0
terminal_n_bp = 0
short_bp = 0

with open(output_fa, "w") as out:

    for header, seq in read_fasta(input_fa):

        input_records += 1
        input_bp += len(seq)

        seqid = header.split()[0]

        # Remove whole-contig contamination
        if seqid in exclude:
            print(f"EXCLUDE\t{seqid}\t{len(seq)} bp")
            excluded_bp += len(seq)
            continue

        # Split at internal adaptor
        if seqid in split_regions:

            start, end = split_regions[seqid]

            left_seq = seq[:start-1]
            right_seq = seq[end:]

            adapter_bp += end - start + 1

            suffix = header[len(seqid):]

            pieces = [
                (f"{seqid}a{suffix}", left_seq),
                (f"{seqid}b{suffix}", right_seq)
            ]

            print(
                f"SPLIT\t{seqid}\t"
                f"remove {start}..{end}"
            )

        else:
            pieces = [(header, seq)]

        for new_header, piece in pieces:

            piece, left_n, right_n = trim_terminal_ns(piece)

            terminal_n_bp += left_n + right_n

            if len(piece) < MIN_LEN:
                print(
                    f"DROP_SHORT\t{new_header.split()[0]}"
                    f"\t{len(piece)} bp"
                )
                short_bp += len(piece)
                continue

            write_record(out, new_header, piece)

            output_records += 1
            output_bp += len(piece)


print()
print("========== SUMMARY ==========")
print(f"Input sequences       : {input_records}")
print(f"Input bp              : {input_bp}")
print(f"Excluded contamination: {excluded_bp} bp")
print(f"Adaptor removed       : {adapter_bp} bp")
print(f"Terminal Ns trimmed   : {terminal_n_bp} bp")
print(f"Short sequence removed: {short_bp} bp")
print(f"Output sequences      : {output_records}")
print(f"Output bp             : {output_bp}")
print("=============================")

PY

echo
echo "Finished:"
ls -lh "$OUTPUT"
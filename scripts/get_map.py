"""
get_map.py

Generate a 4-column space-delimited map file accepted by selscan from ms or VCF input.

Output columns:
1) Chromosome name (default: "chr1")
2) SNP serial number starting at 1
3) Genetic position (coordinate used for integrating iHH decay in iHS)
4) Physical position (bp)

Notes:
- Column 3: genetic position used for iHH integration (here equal to physical position).
- Column 4: physical position used by selscan to determine gap thresholds for EHH truncation.
"""

import gzip

def open_maybe_gz(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")
    

import argparse
import sys

def parse_ms_positions(ms_file, genome_length):
    physical_positions = []
    with open_maybe_gz(ms_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("positions:"):
                parts = line.split()
                rel_pos = parts[1:]  # skip 'positions:'
                physical_positions = [int(float(p) * genome_length) for p in rel_pos]
                break
    if not physical_positions:
        raise ValueError("No positions line found in ms file.")
    return physical_positions

def parse_vcf_positions(vcf_file):
    physical_positions = []
    with open_maybe_gz(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])
            physical_positions.append(pos)
    if not physical_positions:
        raise ValueError("No variant positions found in VCF file.")
    return physical_positions

def write_map(positions, chrom="chr1", out_file="output.map"):
    with open(out_file, "w") as f:
        ## Write header line describing columns (space-delimited)
        # f.write("chrom snp_index genetic_pos physical_pos\n")
        for i, pos in enumerate(positions, start=0):
            f.write(f"{chrom} rs{i} {pos} {pos}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Generate a 4-column (<chrom> <snp_index> <genetic_pos> <physical_pos>) space-delimited map file from ms or VCF input.",
        epilog="Note: Col3=genetic pos for iHH integration (here equals physical pos, assuming constant recombination rate), Col4=physical pos for maximum allowable EHH gap detection."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ms", help="Input ms file with positions line")
    group.add_argument("--vcf", help="Input VCF file")

    parser.add_argument("--len", type=int,
                        help="Genome length for scaling relative positions in ms file (required if --ms is used)")
    parser.add_argument("--chrom", default="chr1", help="Chromosome name to use in output (default: chr1)")
    parser.add_argument("--out", default="output.map", help="Output map file name (default: output.map)")

    args = parser.parse_args()

    if args.ms and args.len is None:
        print("Error: --len is required when using --ms", file=sys.stderr)
        sys.exit(1)

    if args.ms:
        positions = parse_ms_positions(args.ms, args.len)
    else:
        positions = parse_vcf_positions(args.vcf)

    write_map(positions, chrom=args.chrom, out_file=args.out)
    print(f"Map file written to {args.out}")

if __name__ == "__main__":
    main()

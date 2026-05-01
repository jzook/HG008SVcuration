#!/usr/bin/env python3
import argparse
import sys
import re
import pysam
import truvari


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter SVs by size, filter close BND pairs, and decompose INVs into BND pairs."
    )
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")
    parser.add_argument(
        "--min-size", type=int, default=50, help="Minimum SV size in bp (default: 50)"
    )
    return parser.parse_args()


def get_sv_size(record, min_size):
    """Calculates SV size. Returns infinity to bypass filtering if size is unknown."""
    if not record.alts:
        return 0  # Filter out empty ALTs

    alt_str = record.alts[0]
    svtype = record.info.get("SVTYPE", "")

    # Determine if sequences are explicit (not symbolic or breakends)
    is_symbolic = "<" in alt_str or ">" in alt_str or "[" in alt_str or "]" in alt_str

    # 1. Explicit Sequence Logic (Overrides END entirely)
    if not is_symbolic:
        ref_len = len(record.ref)
        alt_len = len(alt_str)

        # Threshold accounts for 1bp padding (e.g., min_size 50 + 1 = 51bp)
        threshold = min_size + 1

        # If BOTH sequences are strictly smaller than the threshold,
        # it's a short variant (SNP/Indel) and should be filtered.
        if ref_len < threshold and alt_len < threshold:
            return 0  # Evaluates to < min_size in main(), triggering the filter

        # Calculate actual size delta based strictly on sequence length
        return abs(ref_len - alt_len)

    # 2. Handle BNDs
    if svtype == "BND" or "[" in alt_str or "]" in alt_str:
        match = re.search(r"[\]\[](.+?):(\d+)[\]\[]", alt_str)
        if match:
            mate_chrom = match.group(1)
            mate_pos = int(match.group(2))

            if mate_chrom == record.chrom:
                return abs(mate_pos - record.pos)
            else:
                return float("inf")  # Keep translocations
        return float("inf")  # Failed to parse BND, do not filter

    # 3. Handle standard SVs using SVLEN
    if "SVLEN" in record.info:
        svlen = record.info.get("SVLEN")
        if isinstance(svlen, (list, tuple)):
            svlen = svlen[0]
        return abs(svlen)

    # 4. Fallback to calculating from END tag (Only for symbolic variants now)
    if "END" in record.info:
        return abs(record.info.get("END") - record.pos)

    # 5. Final Fallback: If it's symbolic and lacks SVLEN/END, DO NOT FILTER.
    return float("inf")


def main():
    args = parse_args()

    try:
        vcf_in = pysam.VariantFile(args.input)
    except Exception as e:
        sys.exit(f"Error opening input file: {e}")

    if "MATEID" not in vcf_in.header.info:
        vcf_in.header.info.add("MATEID", ".", "String", "ID of mate breakend")

    # Save header and chromosome order before processing
    header = vcf_in.header
    chrom_order = {contig: idx for idx, contig in enumerate(header.contigs)}

    stats = {"inv_decomposed": 0, "filtered_small": 0, "passed": 0}

    # Collect all records first, then sort before writing
    output_records = []

    for record in vcf_in:
        # Pass the args.min_size to our updated sizing function
        if get_sv_size(record, args.min_size) < args.min_size:
            stats["filtered_small"] += 1
            continue

        # Decompose INVs
        if "INV" in record.alts[0] or record.info.get("SVTYPE") == "INV":
            stats["inv_decomposed"] += 1
            trec = truvari.VariantRecord(record)

            bnd_records = trec.decompose()

            for bnd in bnd_records:
                # Get the underlying pysam record from Truvari VariantRecord
                if hasattr(bnd, 'get_record'):
                    output_records.append(bnd.get_record())
                else:
                    output_records.append(bnd)
        else:
            stats["passed"] += 1
            output_records.append(record)

    vcf_in.close()

    # Sort records by chromosome and position
    print("Sorting records...")

    def sort_key(rec):
        chrom_idx = chrom_order.get(rec.contig, 999999)
        return (chrom_idx, rec.pos)

    output_records.sort(key=sort_key)

    # Write sorted records
    output_mode = "wz" if args.output.endswith('.gz') else "w"
    vcf_out = pysam.VariantFile(args.output, output_mode, header=header)

    for record in output_records:
        vcf_out.write(record)

    vcf_out.close()

    print("--- Processing Summary ---")
    print(f"SVs filtered (< {args.min_size}bp): {stats['filtered_small']}")
    print(f"INVs decomposed:         {stats['inv_decomposed']}")
    print(f"Other SVs passed:        {stats['passed']}")


if __name__ == "__main__":
    main()

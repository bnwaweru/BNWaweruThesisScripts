#!/usr/bin/env python3

import csv
import re
import argparse
from collections import defaultdict

HEADER = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

QSEQID_RE = re.compile(
    r'^(?P<snp_chr>[^_]+)_(?P<snp_pos>\d+)::(?P<region_chr>[^:]+):(?P<region_start>\d+)-(?P<region_end>\d+)$'
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter BLAST hits for genome-specific KASP SNP candidates."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="BLAST outfmt 6 file without header"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV with SNP classification"
    )
    parser.add_argument(
        "--detail-output", default=None,
        help="Optional detailed per-hit output TSV"
    )
    parser.add_argument(
        "--locus-chr", default=None,
        help="Only keep SNPs from this chromosome/contig in qseqid, e.g. chr5B"
    )
    parser.add_argument(
        "--locus-start", type=int, default=None,
        help="Only keep SNPs with position >= this value"
    )
    parser.add_argument(
        "--locus-end", type=int, default=None,
        help="Only keep SNPs with position <= this value"
    )
    parser.add_argument(
        "--min-offtarget-length", type=int, default=180,
        help="Minimum off-target alignment length to consider strong (default: 180)"
    )
    parser.add_argument(
        "--min-offtarget-pident", type=float, default=90.0,
        help="Minimum off-target percent identity to consider strong (default: 90.0)"
    )
    parser.add_argument(
        "--min-bitscore-ratio", type=float, default=0.70,
        help="Reject if off-target bitscore >= this fraction of best bitscore (default: 0.70)"
    )
    parser.add_argument(
        "--review-bitscore-ratio", type=float, default=0.50,
        help="Review if off-target bitscore >= this fraction of best bitscore (default: 0.50)"
    )
    return parser.parse_args()

def parse_qseqid(qseqid):
    """
    Expected format:
    chr5B_4101::chr5B:4000-4201
    """
    m = QSEQID_RE.match(qseqid)
    if not m:
        return None
    return {
        "snp_chr": m.group("snp_chr"),
        "snp_pos": int(m.group("snp_pos")),
        "region_chr": m.group("region_chr"),
        "region_start": int(m.group("region_start")),
        "region_end": int(m.group("region_end")),
    }

def load_hits(path):
    hits_by_query = defaultdict(list)

    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row_num, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) != 12:
                raise ValueError(
                    f"Line {row_num}: expected 12 columns, found {len(row)}"
                )

            rec = dict(zip(HEADER, row))

            # Convert numeric fields
            rec["pident"] = float(rec["pident"])
            rec["length"] = int(rec["length"])
            rec["mismatch"] = int(rec["mismatch"])
            rec["gapopen"] = int(rec["gapopen"])
            rec["qstart"] = int(rec["qstart"])
            rec["qend"] = int(rec["qend"])
            rec["sstart"] = int(rec["sstart"])
            rec["send"] = int(rec["send"])
            rec["evalue"] = float(rec["evalue"])
            rec["bitscore"] = float(rec["bitscore"])

            parsed = parse_qseqid(rec["qseqid"])
            if parsed is None:
                raise ValueError(
                    f"Could not parse qseqid format on line {row_num}: {rec['qseqid']}"
                )
            rec.update(parsed)

            hits_by_query[rec["qseqid"]].append(rec)

    return hits_by_query

def in_locus(rec, locus_chr=None, locus_start=None, locus_end=None):
    if locus_chr is not None and rec["snp_chr"] != locus_chr:
        return False
    if locus_start is not None and rec["snp_pos"] < locus_start:
        return False
    if locus_end is not None and rec["snp_pos"] > locus_end:
        return False
    return True

def looks_like_self_hit(hit):
    """
    Self-hit definition:
    - subject chromosome matches source interval chromosome
    - aligned subject coordinates overlap the source interval
    """
    if hit["sseqid"] != hit["region_chr"]:
        return False

    smin = min(hit["sstart"], hit["send"])
    smax = max(hit["sstart"], hit["send"])
    rmin = hit["region_start"]
    rmax = hit["region_end"]

    overlap = not (smax < rmin or smin > rmax)
    return overlap

def classify_query(hits, args):
    """
    1. Find best hit by bitscore
    2. Prefer the best self-hit if present
    3. Evaluate all other hits as off-target
    """
    hits_sorted = sorted(hits, key=lambda x: (-x["bitscore"], x["evalue"]))

    self_hits = [h for h in hits_sorted if looks_like_self_hit(h)]
    best_hit = hits_sorted[0]
    best_self_hit = self_hits[0] if self_hits else None

    # Use best self-hit as the reference if present; otherwise use best hit
    anchor_hit = best_self_hit if best_self_hit is not None else best_hit
    best_bitscore = anchor_hit["bitscore"]
    query_len = max(h["qend"] for h in hits_sorted)

    strong_offtargets = []
    moderate_offtargets = []

    for h in hits_sorted:
        if h is anchor_hit:
            continue

        # Any other self-like hit on same interval is ignored only if it is essentially the same location
        same_subject = h["sseqid"] == anchor_hit["sseqid"]
        same_coords = (
            min(h["sstart"], h["send"]) == min(anchor_hit["sstart"], anchor_hit["send"])
            and max(h["sstart"], h["send"]) == max(anchor_hit["sstart"], anchor_hit["send"])
        )
        if same_subject and same_coords:
            continue

        bitscore_ratio = h["bitscore"] / best_bitscore if best_bitscore > 0 else 0.0
        h["bitscore_ratio"] = bitscore_ratio

        is_strong = (
            h["length"] >= args.min_offtarget_length and
            h["pident"] >= args.min_offtarget_pident and
            bitscore_ratio >= args.min_bitscore_ratio
        )

        is_moderate = (
            h["length"] >= max(100, int(0.6 * query_len)) and
            h["pident"] >= 85.0 and
            bitscore_ratio >= args.review_bitscore_ratio
        )

        if is_strong:
            strong_offtargets.append(h)
        elif is_moderate:
            moderate_offtargets.append(h)

    if strong_offtargets:
        decision = "REJECT"
        reason = f"{len(strong_offtargets)} strong off-target hit(s)"
    elif moderate_offtargets:
        decision = "REVIEW"
        reason = f"{len(moderate_offtargets)} moderate off-target hit(s)"
    else:
        decision = "KEEP"
        reason = "No strong off-target hits"

    return {
        "decision": decision,
        "reason": reason,
        "anchor_hit": anchor_hit,
        "best_hit": best_hit,
        "strong_offtargets": strong_offtargets,
        "moderate_offtargets": moderate_offtargets,
    }

def main():
    args = parse_args()
    hits_by_query = load_hits(args.input)

    summary_rows = []
    detail_rows = []

    for qseqid, hits in hits_by_query.items():
        representative = hits[0]

        if not in_locus(
            representative,
            locus_chr=args.locus_chr,
            locus_start=args.locus_start,
            locus_end=args.locus_end
        ):
            continue

        result = classify_query(hits, args)
        anchor = result["anchor_hit"]

        summary_rows.append({
            "qseqid": qseqid,
            "snp_chr": representative["snp_chr"],
            "snp_pos": representative["snp_pos"],
            "region_chr": representative["region_chr"],
            "region_start": representative["region_start"],
            "region_end": representative["region_end"],
            "decision": result["decision"],
            "reason": result["reason"],
            "best_subject": anchor["sseqid"],
            "best_pident": anchor["pident"],
            "best_length": anchor["length"],
            "best_sstart": anchor["sstart"],
            "best_send": anchor["send"],
            "best_bitscore": anchor["bitscore"],
            "n_hits_total": len(hits),
            "n_strong_offtargets": len(result["strong_offtargets"]),
            "n_moderate_offtargets": len(result["moderate_offtargets"]),
        })

        for h in hits:
            bitscore_ratio = h["bitscore"] / anchor["bitscore"] if anchor["bitscore"] > 0 else 0.0
            detail_rows.append({
                "qseqid": qseqid,
                "snp_chr": h["snp_chr"],
                "snp_pos": h["snp_pos"],
                "region_chr": h["region_chr"],
                "region_start": h["region_start"],
                "region_end": h["region_end"],
                "sseqid": h["sseqid"],
                "pident": h["pident"],
                "length": h["length"],
                "qstart": h["qstart"],
                "qend": h["qend"],
                "sstart": h["sstart"],
                "send": h["send"],
                "evalue": h["evalue"],
                "bitscore": h["bitscore"],
                "bitscore_ratio_to_anchor": round(bitscore_ratio, 4),
                "is_anchor_hit": "YES" if h is anchor else "NO",
                "looks_like_self_hit": "YES" if looks_like_self_hit(h) else "NO",
                "final_decision": result["decision"],
            })

    # Write summary
    summary_header = [
        "qseqid", "snp_chr", "snp_pos", "region_chr", "region_start", "region_end",
        "decision", "reason",
        "best_subject", "best_pident", "best_length", "best_sstart", "best_send", "best_bitscore",
        "n_hits_total", "n_strong_offtargets", "n_moderate_offtargets"
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=summary_header, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    # Optional detailed output
    if args.detail_output:
        detail_header = [
            "qseqid", "snp_chr", "snp_pos", "region_chr", "region_start", "region_end",
            "sseqid", "pident", "length", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore",
            "bitscore_ratio_to_anchor", "is_anchor_hit", "looks_like_self_hit", "final_decision"
        ]
        with open(args.detail_output, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=detail_header, delimiter="\t")
            writer.writeheader()
            writer.writerows(detail_rows)

if __name__ == "__main__":
    main()
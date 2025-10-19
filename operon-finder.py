#!/usr/bin/env python3
import sys
import argparse

def parse_line(line):
    """
    Expected columns (whitespace-delimited):
    1: query_id_with_species   e.g., 'WP_023094178.1#1|Pseudomonas_aeruginosa'
    2: (ignored)
    3: (ignored)
    4: strand (+/-)
    5-7: (ignored)
    8: start
    9: stop
    10: neighbor_protein_id    e.g., 'WP_047925767.1#1' OR the query itself for the self row
    11: contig                 e.g., 'NZ_...'
    12: assembly               e.g., 'GCF_...'
    """
    parts = line.strip().split()
    if len(parts) < 12:
        return None  # skip malformed rows
    obj = {
        "q_full": parts[0],                      # col 1
        "q_id": parts[0].split("|")[0],          # query protein id without species
        "q_species": (parts[0].split("|",1)[1] if "|" in parts[0] else ""),
        "strand": parts[3],                      # col 4
        "start": int(parts[7]),                  # col 8
        "stop":  int(parts[8]),                  # col 9
        "prot_id": parts[9],                     # col 10 (neighbor/self)
        "contig": parts[10],                     # col 11
        "assembly": parts[11],                   # col 12
    }
    if obj["start"] > obj["stop"]:
        obj["start"], obj["stop"] = obj["stop"], obj["start"]
    return obj

def blocks_by_query(lines):
    """
    Group consecutive rows into blocks by the value of column 1 (q_full).
    This matches your file where each context block repeats the same col-1 entry,
    then switches to a new one. (No need for blank-line separators.)
    """
    current_key = None
    current = []
    for line in lines:
        if not line.strip():
            # ignore blank lines entirely
            continue
        row = parse_line(line)
        if row is None:
            continue
        key = row["q_full"]
        if current_key is None:
            current_key = key
        if key != current_key:
            yield current_key, current
            current_key = key
            current = []
        current.append(row)
    if current:
        yield current_key, current

def build_operons(rows, max_gap=50, require_same_strand=False):
    """
    Sort genes by start and chain them into operon clusters:
    two adjacent genes belong to the same operon if
    (next.start - prev.stop) <= max_gap (and optionally same strand).
    Returns list of operons; each operon is a list of row dicts.
    """
    rows_sorted = sorted(rows, key=lambda r: (r["contig"], r["start"], r["stop"]))
    operons = []
    current = []
    for r in rows_sorted:
        if not current:
            current = [r]
            continue
        prev = current[-1]
        # break operon if contig changes
        if r["contig"] != prev["contig"]:
            operons.append(current)
            current = [r]
            continue
        intergenic = r["start"] - prev["stop"]
        same_strand_ok = (r["strand"] == prev["strand"]) if require_same_strand else True
        if intergenic <= max_gap and same_strand_ok:
            current.append(r)
        else:
            operons.append(current)
            current = [r]
    if current:
        operons.append(current)
    return operons

def main():
    ap = argparse.ArgumentParser(
        description="Find contexts where the query protein (col 1) is in an operon "
                    "with other proteins based on start/stop (cols 8/9) and IDs (col 10)."
    )
    ap.add_argument("infile", help="Input text file from FlaGs2 context table")
    ap.add_argument("--max-gap", type=int, default=50,
                    help="Maximum intergenic distance (bp) to chain genes into an operon (default: 50)")
    ap.add_argument("--require-same-strand", action="store_true",
                    help="If set, only chain genes on the same strand into operons")
    ap.add_argument("--out", default="-",
                    help="Output TSV file (default: stdout)")
    args = ap.parse_args()

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    # Header
    out.write("\t".join([
        "query_full", "query_id", "query_species", "assembly", "contig",
        "operon_start", "operon_end", "operon_size",
        "members_protein_ids", "members_coords", "query_in_operon"
    ]) + "\n")

    with open(args.infile) as fh:
        for q_full, rows in blocks_by_query(fh):
            if not rows:
                continue
            query_id = rows[0]["q_id"]  # same for the whole block
            assembly = rows[0]["assembly"]
            # build operons
            operons = build_operons(rows, max_gap=args.max_gap,
                                    require_same_strand=args.require_same_strand)

            # find the operon that contains the query gene (self row has prot_id == query_id)
            query_operons = []
            for op in operons:
                prot_ids = {r["prot_id"] for r in op}
                if query_id in prot_ids:
                    query_operons.append(op)

            # If multiple contigs/operons contain a self row (unusual), report all
            for op in query_operons:
                # Only keep if there are other members besides the query itself
                members = [r for r in op]
                # Count non-self members:
                non_self = [r for r in members if r["prot_id"] != query_id]
                if len(non_self) == 0:
                    continue  # query alone; not "associated with other proteins"

                contig = members[0]["contig"]
                op_start = min(r["start"] for r in members)
                op_end   = max(r["stop"]  for r in members)
                members_ids = ",".join(r["prot_id"] for r in members)
                members_coords = ",".join(f"{r['prot_id']}:{r['start']}-{r['stop']}" for r in members)

                out.write("\t".join([
                    q_full,
                    query_id,
                    rows[0]["q_species"],
                    assembly,
                    contig,
                    str(op_start),
                    str(op_end),
                    str(len(members)),
                    members_ids,
                    members_coords,
                    "YES"
                ]) + "\n")

    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()

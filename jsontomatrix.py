#!/usr/bin/env python3

import sys
import json
import csv

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} session.json output.csv")
    sys.exit(1)

session_file = sys.argv[1]
out_file = sys.argv[2]

# Load clinker session JSON
with open(session_file) as f:
    data = json.load(f)

clusters = data["clusters"]   # dict: cluster_id -> {uid, name, loci}
loci = data["loci"]           # dict: locus_id -> {uid, name, genes}
genes = data["genes"]         # dict: gene_id -> {uid, label, names, ...}
groups = data["groups"]       # list of group dicts
cluster_order = data["order"] # list of cluster_ids in display order

# Get cluster (genome) names in a stable order
cluster_names = [clusters[cid]["name"] for cid in cluster_order]

# Map each gene UID to the cluster (genome) it belongs to
gene_to_cluster = {}
for cid in cluster_order:
    cname = clusters[cid]["name"]
    for locus_id in clusters[cid]["loci"]:
        locus = loci[locus_id]
        for gid in locus["genes"]:
            gene_to_cluster[gid] = cname

rows = []

# ---------- 1) Homology groups from clinker ----------
for group in groups:
    group_label = group["label"]      # e.g. "Group 0"
    group_gene_ids = group["genes"]   # list of gene IDs in this group

    # presence/absence per cluster
    presence = {c: 0 for c in cluster_names}
    for gid in group_gene_ids:
        cname = gene_to_cluster.get(gid)
        if cname is not None:
            presence[cname] = 1

    present_clusters = [c for c in cluster_names if presence[c] == 1]

    if len(present_clusters) == len(cluster_names):
        pattern = "all_present"
    elif len(present_clusters) == 1:
        pattern = f"{present_clusters[0]}_only"
    elif len(present_clusters) == 0:
        pattern = "none"
    else:
        pattern = "partial_" + "_".join(present_clusters)

    # representative gene info
    rep_gene = genes[group_gene_ids[0]]
    rep_name = rep_gene["names"].get("locus_tag", rep_gene["label"])
    rep_product = rep_gene["names"].get("product", "")

    row = [
        group_label,
        rep_name,
        rep_product,
        pattern
    ] + [presence[c] for c in cluster_names]

    rows.append(row)

# ---------- 2) Singleton genes (not in any group) ----------
all_gene_ids = set(genes.keys())
group_gene_ids = set(gid for g in groups for gid in g["genes"])
singleton_gene_ids = all_gene_ids - group_gene_ids

for gid in singleton_gene_ids:
    g = genes[gid]
    cname = gene_to_cluster[gid]

    # presence: only this genome has it
    presence = {c: 0 for c in cluster_names}
    presence[cname] = 1

    rep_name = g["names"].get("locus_tag", g["label"])
    rep_product = g["names"].get("product", "")

    group_label = f"Singleton_{rep_name}"
    pattern = f"{cname}_only"

    row = [
        group_label,
        rep_name,
        rep_product,
        pattern
    ] + [presence[c] for c in cluster_names]

    rows.append(row)

# ---------- Write CSV ----------
with open(out_file, "w", newline="") as f:
    writer = csv.writer(f)
    header = ["Group", "Representative_gene", "Product", "Pattern"] + cluster_names
    writer.writerow(header)
    writer.writerows(rows)

print(f"Presence/absence table written to {out_file}")

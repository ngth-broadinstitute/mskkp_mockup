#!/usr/bin/env python3
"""Build static JSON endpoints for the dataset-vs-dataset comparison mockup."""

import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE_ROOT = ROOT.parent / "singlecell_ingest_10" / "processed"
PUBLIC = ROOT / "public"
DATA = PUBLIC / "data"
RANDOM_SEED = 7
MAX_ATLAS_CELLS = 2400
MAX_VALUES_PER_GROUP = 180
TOP_MARKERS_PER_CELLTYPE = 35
MAX_GENE_ENDPOINTS = 260

CURATED_GENES = [
    "CXCL12", "LEPR", "PDGFRA", "VCAM1", "COL1A1", "COL2A1", "ACAN", "SOX9",
    "RUNX2", "SP7", "ALPL", "BGLAP", "IBSP", "SOST", "DMP1", "PECAM1", "VWF",
    "KDR", "LYZ", "S100A8", "S100A9", "MPO", "MS4A1", "CD79A", "CD3D", "NKG7",
    "HBB", "GYPA", "PPBP", "PF4",
]

CELLTYPE_COLORS = {
    "osteoblast": "#4E79A7",
    "osteocyte": "#A0CBE8",
    "chondrocyte": "#59A14F",
    "stromal_mesenchymal": "#F28E2B",
    "endothelial": "#76B7B2",
    "myeloid": "#E15759",
    "b_cell": "#B07AA1",
    "t_cell_or_nk": "#EDC948",
    "erythroid_or_mk": "#9C755F",
    "other": "#BAB0AC",
}


def clean_filename(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def stable_hash(*parts):
    text = "|".join(str(p) for p in (RANDOM_SEED, *parts))
    value = 0
    for char in text:
        value = (value * 131 + ord(char)) % 2**32
    return value


def stable_sample(items, cap, *seed_parts):
    items = list(items)
    if len(items) <= cap:
        return items
    keyed = sorted((stable_hash(*seed_parts, item[0] if isinstance(item, tuple) else item), item) for item in items)
    return [item for _, item in keyed[:cap]]


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, separators=(",", ":"), sort_keys=True)
        handle.write("\n")


def read_sample_metadata(path):
    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_marker_genes(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def dataset_id(dataset_name):
    return clean_filename(dataset_name)


def build_dataset_records():
    records = []
    for dataset_dir in sorted(p for p in SOURCE_ROOT.iterdir() if p.is_dir()):
        meta = read_json(dataset_dir / "dataset_metadata.json")
        samples = read_sample_metadata(dataset_dir / "sample_metadata.tsv.gz")
        markers = read_marker_genes(dataset_dir / "marker_genes.tsv")
        did = dataset_id(meta["datasetName"])
        records.append({
            "id": did,
            "dir": dataset_dir,
            "metadata": meta,
            "samples": samples,
            "markers": markers,
        })
    return records


def write_index(records):
    datasets = []
    celltypes = set()
    for record in records:
        meta = record["metadata"]
        counts = Counter(row["cell_type__kp"] for row in record["samples"])
        celltypes.update(counts)
        datasets.append({
            "id": record["id"],
            "datasetName": meta["datasetName"],
            "geoSeries": meta.get("geoSeries"),
            "species": meta.get("species"),
            "organ": meta.get("organ"),
            "tissue": meta.get("tissue"),
            "assay": meta.get("assay"),
            "summary": meta.get("summary"),
            "n_cells": len(record["samples"]),
            "cell_type_counts": dict(sorted(counts.items())),
            "data": {
                "atlas": f"data/datasets/{record['id']}/atlas.json",
                "summary": f"data/datasets/{record['id']}/summary.json",
            },
        })
    human = next((d["id"] for d in datasets if d["species"] == "Homo sapiens"), datasets[0]["id"])
    mouse = next((d["id"] for d in datasets if d["species"] == "Mus musculus"), datasets[-1]["id"])
    write_json(DATA / "index.json", {
        "generated_from": "singlecell_ingest_10/processed",
        "api_kind": "static-json-github-pages",
        "default_left_dataset": human,
        "default_right_dataset": mouse,
        "datasets": datasets,
        "cell_types": sorted(celltypes),
        "cell_type_colors": CELLTYPE_COLORS,
    })


def write_dataset_endpoints(records):
    for record in records:
        rows = record["samples"]
        by_type = defaultdict(list)
        for row in rows:
            by_type[row["cell_type__kp"]].append(row)
        sampled = []
        for cell_type, items in sorted(by_type.items()):
            quota = max(30, round(MAX_ATLAS_CELLS * len(items) / max(1, len(rows))))
            for row in stable_sample(items, min(len(items), quota), record["id"], cell_type, "atlas"):
                sampled.append({
                    "id": row["ID"],
                    "x": round(float(row["X"]), 4),
                    "y": round(float(row["Y"]), 4),
                    "cell_type": row["cell_type__kp"],
                    "biosample_id": row["biosample_id"],
                })
        if len(sampled) > MAX_ATLAS_CELLS:
            sampled = stable_sample(sampled, MAX_ATLAS_CELLS, record["id"], "atlas-final")
        meta = record["metadata"]
        counts = Counter(row["cell_type__kp"] for row in rows)
        write_json(DATA / "datasets" / record["id"] / "atlas.json", {
            "dataset_id": record["id"],
            "datasetName": meta["datasetName"],
            "n_cells_total": len(rows),
            "n_cells_exported": len(sampled),
            "downsampling": f"Stratified by cell_type__kp to at most {MAX_ATLAS_CELLS} cells.",
            "cells": sampled,
        })
        write_json(DATA / "datasets" / record["id"] / "summary.json", {
            "dataset_id": record["id"],
            "metadata": {k: meta.get(k) for k in ["datasetName", "geoSeries", "species", "organ", "tissue", "assay", "summary", "pmid", "doi"]},
            "n_cells": len(rows),
            "cell_type_counts": dict(sorted(counts.items())),
        })


def collect_gene_panel(records):
    genes = []
    seen = set()
    for gene in CURATED_GENES:
        seen.add(gene.upper())
        genes.append(gene)
    for record in records:
        for row in record["markers"]:
            if int(float(row.get("rank") or 9999)) > TOP_MARKERS_PER_CELLTYPE:
                continue
            gene = row["gene"]
            key = gene.upper()
            if key not in seen:
                seen.add(key)
                genes.append(gene)
            if len(genes) >= MAX_GENE_ENDPOINTS:
                return genes
    return genes


def read_expression_rows(norm_counts_path, wanted):
    wanted_upper = {gene.upper(): gene for gene in wanted}
    out = {}
    with gzip.open(norm_counts_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene = fields[0]
            key = gene.upper()
            if key not in wanted_upper:
                continue
            out[wanted_upper[key]] = [float(v) if v else 0.0 for v in fields[1:]]
    return cells, out


def values_by_celltype(samples, values, dataset, gene, cell_type):
    pairs = [(samples[i]["ID"], values[i]) for i in range(min(len(samples), len(values))) if samples[i]["cell_type__kp"] == cell_type]
    sampled = stable_sample(pairs, MAX_VALUES_PER_GROUP, dataset, gene, cell_type)
    return [round(float(value), 4) for _, value in sampled]


def summarize_values(values):
    if not values:
        return {"n": 0, "avg_expression": 0, "pct_expressing": 0, "median": 0}
    ordered = sorted(values)
    mid = len(ordered) // 2
    median = ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2
    return {
        "n": len(values),
        "avg_expression": round(sum(values) / len(values), 5),
        "pct_expressing": round(sum(1 for value in values if value > 0) / len(values), 5),
        "median": round(median, 5),
    }


def write_gene_and_celltype_endpoints(records):
    gene_panel = collect_gene_panel(records)
    gene_payloads = {
        gene: {
            "gene": gene,
            "downsampling_note": f"Violin values sampled up to {MAX_VALUES_PER_GROUP} cells per dataset per cell type with fixed seed {RANDOM_SEED}.",
            "datasets": {},
        }
        for gene in gene_panel
    }
    celltype_gene_candidates = defaultdict(set)
    celltype_payloads = defaultdict(lambda: {"cell_type": None, "genes": {}})

    for record in records:
        did = record["id"]
        samples = record["samples"]
        cell_types = sorted(set(row["cell_type__kp"] for row in samples))
        cells, expr = read_expression_rows(record["dir"] / "norm_counts.tsv.gz", gene_panel)
        if cells != [row["ID"] for row in samples]:
            raise RuntimeError(f"Cell order mismatch for {did}")
        for row in record["markers"]:
            if int(float(row.get("rank") or 9999)) <= 12:
                celltype_gene_candidates[row["cell_type__kp"]].add(row["gene"])
        for gene, values in expr.items():
            dataset_rows = []
            for cell_type in cell_types:
                full_values = [values[i] for i, sample in enumerate(samples) if sample["cell_type__kp"] == cell_type]
                dataset_rows.append({
                    "cell_type": cell_type,
                    "summary": summarize_values(full_values),
                    "values": values_by_celltype(samples, values, did, gene, cell_type),
                })
            gene_payloads[gene]["datasets"][did] = dataset_rows

    for gene, payload in gene_payloads.items():
        write_json(DATA / "genes" / f"{clean_filename(gene)}.json", payload)
    write_json(DATA / "genes" / "index.json", {
        "items": [{"gene": gene, "file": f"data/genes/{clean_filename(gene)}.json"} for gene in gene_panel],
    })

    # Cell-type endpoints reuse gene expression endpoints in the frontend, but
    # keep a compact per-cell-type gene list for pair-specific plotting.
    panel_by_upper = {gene.upper(): gene for gene in gene_panel}
    curated_upper = {gene.upper() for gene in CURATED_GENES}
    for cell_type, genes in sorted(celltype_gene_candidates.items()):
        available = {panel_by_upper[gene.upper()] for gene in genes if gene.upper() in panel_by_upper}
        selected = sorted(available, key=lambda gene: (gene.upper() not in curated_upper, gene.upper()))[:14]
        write_json(DATA / "celltypes" / f"{clean_filename(cell_type)}.json", {
            "cell_type": cell_type,
            "genes": selected,
            "note": "Top marker genes across the ingested datasets for this broad cell label.",
        })


def main():
    records = build_dataset_records()
    write_index(records)
    write_dataset_endpoints(records)
    write_gene_and_celltype_endpoints(records)
    write_json(DATA / "build_info.json", {
        "source": str(SOURCE_ROOT),
        "datasets": [record["id"] for record in records],
        "max_atlas_cells": MAX_ATLAS_CELLS,
        "max_values_per_group": MAX_VALUES_PER_GROUP,
    })


if __name__ == "__main__":
    main()

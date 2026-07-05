# MSKKP Dataset Comparison Mockup

Static GitHub Pages mockup for comparing standardized MSKKP single-cell datasets.

The site is fully static:

- Frontend: `public/index.html`, `public/styles.css`, `public/app.js`
- Static API endpoints: `public/data/**/*.json`
- Data builder: `build_static_api.py`

The current default comparison is human `GSE196678` vs mouse `GSE122465`, with `GSE106236` also available in the dataset selectors.

Dataset selectors use descriptive study labels in the UI while retaining the canonical dataset IDs in the JSON endpoints.

Cell-type comparison endpoints contain one summary row per retained gene for each broad cell type. The scatter plot uses these summary-only values:

- x-axis: average log-normalized expression in the left dataset
- y-axis: average log-normalized expression in the right dataset
- gene matching: case-insensitive gene symbol key

Per-cell sampled distributions are still generated only for the gene-comparison violin endpoints. Those values are sampled up to 180 cells per dataset per cell type with fixed seed `7`; groups below the cap keep all cells.

## Rebuild Data

From the parent project directory:

```bash
python3 dataset_compare_pages/build_static_api.py
```

The builder reads canonical ingested outputs from:

`../singlecell_ingest_10/processed`

## Deployment

GitHub Pages is configured to serve the `gh-pages` branch. After rebuilding `public/`, copy the static files to that branch and push it.

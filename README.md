# MSKKP Dataset Comparison Mockup

Static GitHub Pages mockup for comparing standardized MSKKP single-cell datasets.

The site is fully static:

- Frontend: `public/index.html`, `public/styles.css`, `public/app.js`
- Static API endpoints: `public/data/**/*.json`
- Data builder: `build_static_api.py`

The current default comparison is human `GSE196678` vs mouse `GSE122465`, with `GSE106236` also available in the dataset selectors.

## Rebuild Data

From the parent project directory:

```bash
python3 dataset_compare_pages/build_static_api.py
```

The builder reads canonical ingested outputs from:

`../singlecell_ingest_10/processed`

## Deployment

This repository includes a GitHub Actions workflow in `.github/workflows/pages.yml` that uploads `public/` to GitHub Pages.

After pushing to GitHub, enable Pages with GitHub Actions as the source if it is not already enabled.

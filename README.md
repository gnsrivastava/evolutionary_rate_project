# evolutionary_rate_project
Evolutionary rate calculator based on published molecular clock data

## Input support

The pipeline accepts TSV/CSV manifests with flexible column aliases. For mixed-source manifests, it now preserves and uses:
- `species_query` / `species_name`
- `assembly_accession`
- `ncbi_taxon_id`
- `mapping_level`
- `mapping_description`
- `acdb_strain_designation` / `ncbi_strain`
- `selection_method`

## BV-BRC and non-BV-BRC behavior

- If BV-BRC genome IDs are resolved, the pipeline enriches results with BV-BRC metadata, AMR phenotypes, and specialty genes.
- If a record is not found in BV-BRC, it is still processed with taxonomy-driven rate lookup (static/dynamic/default), and BV-BRC-only fields are left as missing values.
- Output includes resolution diagnostics:
  - `resolution_source`
  - `bvbrc_available`
  - `proxy_used`
  - `skip_reason`

## Strict mode

Use `--strict-bvbrc` to keep only records that have BV-BRC metadata support.

## Composite score behavior

- `composite_score` is computed with fixed absolute scaling and is independent of how many species are in the same run.
- `rate_per_site_per_year` is log-scaled using fixed bounds (`1e-9` to `1e-4`).
- `amr_genes` is scaled using a fixed cap (`20` genes).
- `resistance_fraction` is clipped to `[0, 1]` (missing values use `0.5`).
- Use `--normalization-mode absolute` (default) for fixed absolute scaling.
- Use `--normalization-mode batch` to compute dataset-relative (batch-normalized) component scores before `composite_score` aggregation.

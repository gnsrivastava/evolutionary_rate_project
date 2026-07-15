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

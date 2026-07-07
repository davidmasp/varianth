# Running `varianth gene2protein` with MANE (GRCh38 v1.4)

This note describes a practical workflow to generate and post-process the full nonsynonymous mutation table using the `gene2protein` subcommand and [MANE](https://ncbiinsights.ncbi.nlm.nih.gov/2024/10/28/mane-v1-4-mane-select-non-coding-genes/) inputs.

## Goal

Create a complete table of possible nonsynonymous SNVs for MANE proteins, then produce two sorted views:

- genomic-position sorted (best for locus/chromosome workflows)
- protein-position sorted (best for protein/codon workflows)

## Inputs

- `MANE.GRCh38.v1.4.ensembl_genomic.gff.gz`
- `genome.fa` (+ `genome.fa.fai`) - for example from iGenomes
- `MANE.GRCh38.v1.4.ensembl_protein.faa` (+ `.fai`) - ⏰ Note that this file needs to be de-compressed and indexed by `samtools faidx` before procedding. 

## 1) Run `gene2protein`

```bash
# 1 min
varianth gene2protein \
  --gff-path MANE.GRCh38.v1.4.ensembl_genomic.gff.gz \
  --genome-fasta-path genome.fa \
  --proteome-fasta-path MANE.GRCh38.v1.4.ensembl_protein.faa \
  --output-prefix tables/all_mutations
```

This generates:

- `tables/all_mutations.tsv`
- `tables/all_mutations.json`

## 2) Sort output for different access patterns

Set locale for byte-wise stable sorting and better performance:

```bash
export LC_ALL=C
```

Sort by genomic coordinates (`chrom`, `genome_position`):

```bash
# ~ 1min
sort --parallel=4 -S 8G -k1,1 -k2,2n tables/all_mutations.tsv > tables/all_mutations_gsorted.tsv
```

Sort by protein coordinates (`protein_id`, `protein_position`):

```bash
# ~ 1min
sort --parallel=4 -S 8G -k5,5 -k7,7n tables/all_mutations.tsv > tables/all_mutations_psorted.tsv
```

## 3) Compression and indexing

Compress and index the genomic-sorted file using bgzip:

```bash
bgzip tables/all_mutations_gsorted.tsv
tabix -s 1 -b 2 -e 2 all_mutations_gsorted.tsv.gz
```

This allows for random access at any genome coordinate (as in `tabix all_mutations_gsorted.tsv.gz chr2:47482936-47482976`).

We can do the same with protein coordinates.

```bash
bgzip all_mutations_psorted.tsv
tabix -s 5 -b 7 -e 7 all_mutations_psorted.tsv.gz
```

This also allows for random access at any protein context (as in `tabix all_mutations_psorted.tsv.gz ENSP00000233146.2:931-932`).


## Why each output file is useful

- `all_mutations.tsv`: Canonical raw output from `gene2protein`; best as a reproducible source table.
- `all_mutations.json`: Run diagnostics/QA (`total_proteins`, `successful_count`, `failed_count`, failed IDs/errors, elapsed time).
- `all_mutations_gsorted.tsv.gz` + index: Best for coordinate-based operations (chromosome/position joins, region filters, genomic overlap workflows).
- `all_mutations_psorted.tsv` + index: Best for protein-centric analyses (per-protein scans, amino-acid position joins, codon-level aggregations).

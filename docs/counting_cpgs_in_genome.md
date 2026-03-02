# Counting CpGs in the Genome with `varianth kcount`

This note describes one practical way to calculate the number of CpG dinucleotides in a reference genome using 3-mer counts.

## Goal

Count all 3-mers (`k=3`) in the genome, then sum contexts of the form `NCG` (where `N` is any base).  
Each `NCG` contains one `CG` at positions 2-3. In this workflow we multiply by 2 because both strands are considered.

## 1) Generate k-mer counts

```bash
varianth kcount -K 3 genome.fa --skip-ambiguous > kmer_genome.json
```

Why this file is useful:

- `kmer_genome.json`: reusable per-contig 3-mer count table for CpG estimation and other context analyses.

## 2) Compute CpG calculation in R

```r
library(jsonlite)
library(magrittr)
library(ggplot2)

path = "kmer_genome.json"

dat = jsonlite::fromJSON(path)
dat_df = dat[["counts"]] %>%
  purrr::map2_df(names(.), function(x, name){
    data.frame(
      kmer = names(x),
      count = as.integer(unlist(x)),
      genome = name,
      stringsAsFactors = FALSE
    )
  })

ctx_in = helperMut::make_set(x = "NCG>T") %>%
  stringr::str_extract("[:alnum:]+(?=>T)")

ncg_df = dat_df %>%
  dplyr::group_by(kmer) %>%
  dplyr::summarise(total_count = sum(count)) %>%
  dplyr::filter(
    kmer %in% ctx_in
  )

## Multiply by 2 to account for both strands.
scales::comma(2 * sum(ncg_df$total_count))
```

Why these objects are useful:

- `dat_df`: tidy long-format k-mer table suitable for plotting and downstream filtering.
- `ctx_in`: the specific context set used for `NCG` extraction.
- `ncg_df`: aggregated counts for CpG-containing trinucleotide contexts.
- final scalar (`2 * sum(...)`): genome-wide CpG estimate.

## Notes

- `--skip-ambiguous` avoids failures from ambiguous bases (e.g. `N`) and excludes ambiguous k-mers.
- Keep `K=3` for this specific approach; larger `K` values are useful for broader context profiling but are not needed for direct CpG counting.

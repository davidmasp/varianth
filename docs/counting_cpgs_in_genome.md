



```bash
varianth kcount -K 3 genome.fa --skip-ambiguous > kmer_genome.json
```


```R
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

## i need to multiply because every NCG has to CG dinucleotides.
scales::comma(2*sum(ncg_df$total_count))

```


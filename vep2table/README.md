

```R

library(magrittr)
library(readr)

# params ------------------------------------------------------------------

path = "/path/to/test.csv"

col_types = cols(
  CLIN_SIG = col_character(),
  MANE_PLUS_CLINICAL = col_character(),
  PUBMED = col_character(),
  CDS_position = col_character(),
  Protein_position = col_character(),
  miRNA = col_character(),
  HGVS_OFFSET = col_character()
)

dat = readr::read_delim(path, delim = "|", col_names = T,
                        col_types = col_types)

dat$Consequence %>% table()
```




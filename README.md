# TransAT
TransAT is a package within the R programming language for converting RefSeq ID to Ensembl ID and mapping the genomic position, followed by providing variant allele frequencies and gene annotations.

**Install the package**
```{r setup}
devtools::install_github("ShihChingYu/TransAT", force=T)
library(TransAT)
```

## Usage of convet_tanscriptID()

**Load the data for convet_tanscriptID()**

- load .csv of refseq ID
```{r}
dat<-read.csv(system.file("extdata", "convertID_refseq_data.csv", package = "TransAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
dat
```

- load .csv of UCSC ID
```{r}
dat_ucsc<-read.csv(system.file("extdata", "convertID_ucsc_data.csv", package = "TransAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
dat_ucsc
```

- load .bed of refseq ID
```{r}
dat_b<-read.table(system.file("extdata", "convertID_refseq_data2.bed", package = "TransAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = "\t", header=T)
dat_b
```

**Excute in convert_transcriptID()**

- convert transcript refseq ID
```{r}
new_dat<-convert_transcriptID(dat, db, dat_filter = "refseq_mrna")
```

- convert transcript ucsc ID
```{r}
new_dat_ucsc<-convert_transcriptID(dat_ucsc, db, dat_filter = "ucsc")
```

**Result from the convert_transcriptID()**
```{r}
data(convertID_result, package = "TransAT")
convertID_result
```

## Usage of pop_freq()

**Load the data for getting population allele frequency**
```{r}
anno_freq_data<-read.csv(system.file("extdata",
                          "anno_freq_data.csv",
                          package = "TransAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
anno_freq_data
```

**Excute in pop_freq()**
```{r}
pop_dat<-pop_freq(anno_freq_data, pop="db_gnomAD_exome_freq")
```

**Result from the pop_freq()**
```{r}
data(pop_result, package = "TransAT")
pop_result
```

This is the [PDF](https://github.com/ShihChingYu/TransAT/blob/master/data-raw/plots.pdf) depicting MAF of each global population with barplots.

# TransAT
TransAT is a package within the R programming language for converting RefSeq ID to Ensembl ID and mapping the genomic position, followed by providing variant allele frequencies and gene annotations.

#Install the package
```{r setup}
devtools::install_github("ShihChingYu/TransAT", force=T)
library(TransAT)
```

#Load the data for convet tanscript ID
```{r}
dat<-read.csv(system.file("extdata", "convertID_refseq_data.csv", package = "TransAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
dat
```

#Excute in convert_transcriptID()
```{r}
new_dat<-convert_transcriptID(dat, db, dat_filter = "refseq_mrna")
```

#Result from the convert_transcriptID()
```{r}
data(convertID_result, package = "TransAT")
convertID_result
```

#Load the data for getting population allele frequency
```{r}
anno_freq_data<-read.csv(system.file("extdata",
                          "anno_freq_data.csv",
                          package = "TransAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
anno_freq_data
```

#Excute in pop_freq()
```{r}
pop_dat<-pop_freq(anno_freq_data, pop="db_gnomAD_exome_freq")
```

#Result from the pop_freq()
```{r}
data(pop_result, package = "TransAT")
pop_result
```

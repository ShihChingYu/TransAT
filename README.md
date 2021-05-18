---
output:
  pdf_document: default
  html_document: default
---
messenger RNA annotation tool (MRAT), that provides (i) SNP ID and genomic position for a user-provided transcript ID from patients, and (ii) allele frequencies for the SNPs from publicly available global populations, through two simple command lines. Two functions, “convert_transcriptID” and “pop_freq”, conduct all the steps, and output intermediate genomic information and final allelic annotation, respectively. It accesses allele frequencies for different global populations from the publicly available 1000 genomes, gnomAD and Taiwan Biobank databases.


### Installation

```{r package}
devtools::install_github("ShihChingYu/MRAT")
library(MRAT)
```

### Note
Please install the package on the Mac or remove the old version of R if you are unable to install it. 

---
title: "README.md"
author: "ChingYuShih"
date: "2021/2/17"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## MRAT package

MRAT is a package within the R programming language for converting RefSeq ID to Ensembl transcripts and mapping the genomic position, followed by providing variant allele frequencies and gene annotations.

To run and install the MRAT package
```{r MRAT}
devtools::install_github("ShihChingYu/MRAT")
library(MRAT)
```


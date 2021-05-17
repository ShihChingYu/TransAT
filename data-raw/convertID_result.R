require(devtools)
library(EnsDb.Hsapiens.v75)
library(ensembldb)
library(biomaRt)
library(IRanges)

db=EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
dat<-read.csv(system.file("extdata",
                          "convertID_refseq_data.csv",
                          package = "MRAT"),
              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")


dat<-as.data.frame(dat)

#preprocess data
str<-as.array(dat$Nucleotide_changes)
dat$CDS_start_loc<-as.integer(unlist(stringr::str_extract_all(str, "[0-9]+")))
dat$ref<-substr(dat$Nucleotide_changes, 1, 1)
dat$alt<-stringr::str_sub(dat$Nucleotide_changes, -1, -1)

#converting NCBI RefSeq transcript ID to Ensembl transcript ID
ensembl <- biomaRt::useEnsembl(biomart = biomart_ens, dataset = dat_ens, GRCh=37)
dat_ensmbl_id<-biomaRt::getBM(attributes=BM_att_ens, filters = dat_filter,
                              values = dat[1], mart= ensembl, uniqueRows = TRUE, useCache = FALSE)
matches <- grep("^NM", dat_ensmbl_id$refseq_mrna, ignore.case = T)
dat_ensmbl_id2<-dat_ensmbl_id[matches,]
dat2<-merge(dat, dat_ensmbl_id2, by.x = colnames(dat[1]), by.y = dat_filter, all.y=T)


#loading the db
#mapping to the genome of the Ensembl transcript
db=EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
dat3<-dat2[!is.na(dat2$ensembl_transcript_id), ]
loc<-dat3$CDS_start_loc
tx_name<-dat3$ensembl_transcript_id
cds <- IRanges::IRanges(start = loc, width = 1, name = tx_name)
cds_tx <- ensembldb::cdsToTranscript(cds, db)
cds_tx_gn<-ensembldb::transcriptToGenome(cds_tx, db)
cds_tx_gn_df<-as.data.frame(cds_tx_gn)
cds_tx_gn_df<-cds_tx_gn_df[order(cds_tx_gn_df$tx_id, cds_tx_gn_df$tx_start),]
com_col<-intersect(colnames(dat3), colnames(cds_tx_gn_df))
cds_tx_gn_df<-cds_tx_gn_df[order(cds_tx_gn_df$tx_id, cds_tx_gn_df$tx_start),]
dat3<-dat3[order(dat3$ensembl_transcript_id, dat3$CDS_start_loc),]
dat4<-cbind(dat3, cds_tx_gn_df)

#retrieve reID from biomart = "snps" and dataset = "hsapiens_snp" based on genomic positions in previous step
SNP_M <- data.frame(CHR = dat4$seqnames, START = dat4$start, END = dat4$end)
coords <- apply(SNP_M, 1, paste, collapse = ":")
snp = biomaRt::useEnsembl(biomart = "snps", dataset = "hsapiens_snp", GRCh=37)
rsid<-biomaRt::getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end'), filters = "chromosomal_region", values = coords, mart = snp, useCache = FALSE, uniqueRows = TRUE)
matches2 <- grep("^rs", rsid$refsnp_id, ignore.case = T)
rsid_onlyrs<-rsid[matches2,]
dat4_rsid<-merge(dat4, rsid_onlyrs, by.x = c("seqnames", "start", "end"), by.y = c("chr_name", "chrom_start", "chrom_end"))
colnames(dat4_rsid)[1]<-"chr"
convertID_result<-dat4_rsid[,c( colnames(dat)[1], "Nucleotide_changes", "CDS_start_loc","ref","alt","chr",
                     "ensembl_transcript_id", "hgnc_symbol", "start", "end", "width", "strand", "exon_id",
                     "exon_rank", "tx_start", "tx_end", "refsnp_id")]

write.csv(convertID_result, "data-raw/convertID_result.csv", row.names = FALSE, quote = FALSE)
usethis::use_data(convertID_result, overwrite = TRUE)

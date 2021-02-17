#' @title Convert transcript IDs between different databases
#' @description Preprocessing, converting and mapping
#' @name convert_transcriptID
#'
#' @param dat a dataframe including Transcript_version, nucleotide "C1768G"
#' refers to ref/CDS position/alt.
#' @param db `EnsDb` object. Default is EnsDb.Hsapiens.v75.
#' @param biomart selection of BioMart database. Default is "ensembl".
#' @param dataset BioMart databases includes many datasets. Choose dataset in the database.
#' Default is "hsapiens_gene_ensembl".
#' @param getBM_attributes defines the values of interests.
#' Default shows "refseq_mrna", "ensembl_transcript_id", "hgnc_symbol".
#' The listAttributes function displays all available attributes in the selected dataset.
#' @param filters a vector of filters to query. Default shows "refseq_mrna".
#' @return a new dataset with converting information
#' @export
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' db=EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
#' dat<-read.csv(system.file("extdata",
#'                           "convertID_data.csv",
#'                           package = "MRAT"),
#'               stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' new_dat<-convert_transcriptID(dat, db, biomart, dataset, getBM_attributes, filters)
#'
convert_transcriptID <- function(dat, db, biomart="ensembl",
                                 dataset="hsapiens_gene_ensembl",
                                 getBM_attributes=c("refseq_mrna", "ensembl_transcript_id", "hgnc_symbol"),
                                 filters="refseq_mrna"){

  #preprocess data
  str<-as.array(dat$Nucleotide_changes)
  dat$CDS_start_loc<-as.integer(unlist(stringr::str_extract_all(str, "[0-9]+")))
  dat$ref<-substr(dat$Nucleotide_changes, 1, 1)
  dat$alt<-stringr::str_sub(dat$Nucleotide_changes, -1, -1)
  dat<-dat[order(dat$Transcript_version, dat$CDS_start_loc),]

  #converting NCBI RefSeq transcript ID to Ensembl transcript ID
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh=37)
  dat_ensmbl_id<-biomaRt::getBM(attributes=c("refseq_mrna", "ensembl_transcript_id", "hgnc_symbol"), filters = "refseq_mrna",
                                values = dat$Transcript_version, mart= ensembl, uniqueRows = TRUE, useCache = FALSE)
  dat_ensmbl_id<-dat_ensmbl_id[order(dat_ensmbl_id$refseq_mrna),]
  dat2<-merge(dat, dat_ensmbl_id, by.x = "Transcript_version", by.y = "refseq_mrna", all=T)

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
  dat3<-dat3[order(dat3$ensembl_transcript_id, dat3$CDS_start_loc),]
  dat4<-cbind(dat3, cds_tx_gn_df)
  dat4
}

#' @title search allel frequency among populations and gene annotation information
#' @description using genomic postion and ref/alt information to query gene annotation and allele frequency
#' @name pop_freq
#'
#' @param dat_ori a dataframe including chromosome, start and end postion, ref/alt nucleotide
#' @param pop selection of table. Default is db_TWB_GWG_freq
#' @return a new dataset with converting information
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata",
#'                           "anno_freq_data.csv",
#'                           package = "MRAT"),
#'               stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' pop_dat<-pop_freq(dat, pop)
#'

pop_freq<-function(ori_data, pop="db_gnomAD_exome_freq"){
  colnames(ori_data)<-c("Chr", "Start", "End", "Ref", "Alt")
  setwd("/work1782/cyshih/package/maf_db")
  sqlite <- RSQLite::dbDriver("SQLite")
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),"maf_db.sqlite3")
  data <- RSQLite::dbGetQuery(conn, 'SELECT * FROM db_gnomAD_exome_freq')

  pop_result <- ori_data %>%
    dplyr::left_join(data , by=c("Chr", "Start", "End", "Ref", "Alt"))

  return(pop_result)
  dbDisconnect(conn)

}

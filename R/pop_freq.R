#' @title search allele frequency among populations and gene annotation information
#' @description using genomic postion and ref/alt information to query gene annotation and allele frequency
#' @name pop_freq
#'
#' @param dat_ori a dataframe including chromosome, start and end postion, ref/alt nucleotide
#' @param pop selection of table. Default is db_TWB_NGS_freq
#' @importFrom magrittr `%>%`
#' @return a new dataset with allele frequency and annotation information
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata",
#'                           "anno_freq_data.csv",
#'                           package = "MRAT"),
#'               stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' pop_dat<-pop_freq(dat, pop="db_TWB_NGS_freq")
#'

pop_freq<-function(dat_ori, pop="db_TWB_NGS_freq"){
  colnames(dat_ori)<-c("Chr", "Start", "End", "Ref", "Alt")
  con = DBI::dbConnect(RMySQL::MySQL(), user='cyshih30', password='bioinfo13579', host = "localhost", dbname="allele_freq_db")

  res <- DBI::dbSendQuery(con, statement=paste("SELECT * FROM", pop))
  data<-DBI::dbFetch(res, n = -1)

  pop_result <- dat_ori %>%
    dplyr::left_join(data , by=c("Chr", "Start", "End", "Ref", "Alt"))

  return(pop_result)
  DBI::dbDisconnect(con)
}




#' @title search allel frequency among populations and gene annotation information
#' @description using genomic postion and ref/alt information to query gene annotation and allele frequency
#' @name pop_freq
#'
#' @param dat_ori a dataframe including chromosome, start and end postion, ref/alt nucleotide
#' @param pop selection of table. Default is db_TWB_GWG_freq
#' @importFrom magrittr `%>%`
#' @return a new dataset with converting information
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata",
#'                           "anno_freq_data.csv",
#'                           package = "MRAT"),
#'               stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' db_file = system.file("extdata", "maf_db.db", package = "MRAT")
#' pop_dat<-pop_freq(dat, pop)
#'

pop_freq<-function(dat_ori, pop="db_gnomAD_exome_freq"){
  colnames(dat_ori)<-c("Chr", "Start", "End", "Ref", "Alt")
  db_file = system.file("extdata", "maf_db.db", package = "MRAT")
  con = DBI::dbConnect(RSQLite::SQLite(), dbname = db_file)
  sqlite <- RSQLite::dbDriver("SQLite")
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),"maf_db.sqlite3")
  data <- RSQLite::dbGetQuery(conn, 'SELECT * FROM db_gnomAD_exome_freq')

  pop_result <- dat_ori %>%
    dplyr::left_join(data , by=c("Chr", "Start", "End", "Ref", "Alt"))

  return(pop_result)
  RSQLite::dbDisconnect(conn)

}

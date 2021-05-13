#' @title search allele frequency among populations and gene annotation information
#' @description using genomic postion and ref/alt information to query gene annotation and allele frequency
#' @name pop_freq
#'
#' @param dat_ori a dataframe including chromosome, start and end postion, ref/alt nucleotide
#' @param pop selection of table. Default is db_gnomAD_exome_freq
#' @importFrom magrittr `%>%`
#' @return a new dataset with allele frequency and annotation information
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata",
#'                           "anno_freq_data.csv",
#'                           package = "MRAT"),
#'               stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' pop_dat<-pop_freq(dat, pop="db_gnomAD_exome_freq")
#'

pop_freq<-function(dat_ori, pop="db_gnomAD_exome_freq"){
  colnames(dat_ori)<-c("Chr", "Start", "End", "Ref", "Alt")
  con = DBI::dbConnect(RMySQL::MySQL(), user='vistor', password='vistor', host = "mrat.cgm.ntu.edu.tw", dbname="allele_freq_db")

  gene_list<-c()
  for (j in 1:nrow(dat_ori)) {
    chr <- dat_ori[j, 1]
    start<- dat_ori[j, 2]
    end<- dat_ori[j, 3]
    ref<-  dat_ori[j, 4]
    alt<-  dat_ori[j, 5]
    query <- paste0("SELECT * FROM `", pop, "` where Chr= '", chr, "' AND
                  Start= '", start, "' AND End= '", end, "' AND Ref= '", ref, "' AND Alt= '", alt, "' ;")
    tmp <- DBI::dbGetQuery(con, query)
    gene_list[[j]]<-tmp
  }
  gene_table<-do.call("rbind", gene_list)

  pop_result <- dat_ori %>%
    dplyr::left_join(gene_table , by=c("Chr", "Start", "End", "Ref", "Alt"))

  return(pop_result)

  #make barplot for each sub-population in the dable
  match_col<-grep("freq$", colnames(pop_result), ignore.case = T)
  pop_dat_forplot<-pop_result[,match_col]

  pdf("allplots.pdf",onefile = TRUE)
  for(i in 1:nrow(pop_dat_forplot)){
    dat<-pop_dat_forplot[i, ]
    dat2<-dat %>% gather(population, frequency, 1:ncol(pop_dat_forplot))
    ggpl<-ggplot2::ggplot(dat2, aes(x=population, y=frequency, color=population)) +
      geom_bar(stat = "identity", fill="white") +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    print(ggpl)}
  dev.off()

  DBI::dbDisconnect(con)
}


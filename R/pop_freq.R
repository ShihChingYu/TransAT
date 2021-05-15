#' @title search allele frequency among populations and gene annotation information
#' @description using genomic postion and ref/alt information to query gene annotation and allele frequency
#' @name pop_freq
#'
#' @param dat_ori a dataframe including chromosome, start and end postion, ref/alt nucleotide
#' @param pop selection of table. Default is db_gnomAD_exome_freq
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @importFrom grDevices pdf dev.off
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

  variants_maf<-list()
  maf<-list()
  for (i in 1:nrow(pop_dat_forplot)){
    variants<-pop_dat_forplot[i, ]
    for (j in 1:ncol(pop_dat_forplot)){
      if(j %% 2 == 1){
        maf[j]<-min(variants[,j], variants[,j+1])
      }
      maf_table<-do.call("cbind", maf)
    }
    variants_maf[[i]]<-maf_table
  }
  var_maf_table<-do.call("rbind", variants_maf)
  colnames(var_maf_table)<-unlist(unique(strsplit(colnames(pop_dat_forplot), "_([^_]*_[^_]*)$")))

  pdf("./barplots.pdf",onefile = TRUE)
  for(i in 1:nrow(var_maf_table)){
    dat<-var_maf_table[i, ]
    dat2<-data.frame(dat)
    ggpl<-ggplot2::ggplot(dat2, aes(x=rownames(dat2), y=dat, color=rownames(dat2))) +
      ggplot2::geom_bar(stat = "identity", fill="white") +xlab('Population') +
      ylab('Frequency') + coord_flip() + theme(legend.position = "none")
    print(ggpl)}
  dev.off()

  DBI::dbDisconnect(con)
}


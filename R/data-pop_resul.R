#' search allele frequency among populations and gene annotation information
#'
#' using genomic postion and ref/alt information to query gene annotation and allele frequency
#'
#' @format A data frame with 3 rows and 28 variables:
#' \describe{
#'   \item{\code{Chr}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{Start}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{End}}{double COLUMN_DESCRIPTION}
#'   \item{\code{Ref}}{character COLUMN_DESCRIPTION}
#'   \item{\code{Alt}}{character COLUMN_DESCRIPTION}
#'   \item{\code{Func.knownGene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{Gene.knownGene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{GeneDetail.knownGene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{ExonicFunc.knownGene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{AAChange.knownGene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_ALL_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_ALL_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_AFR_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_AFR_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_AMR_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_AMR_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_ASJ_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_ASJ_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_EAS_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_EAS_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_FIN_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_FIN_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_NFE_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_NFE_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_SAS_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_SAS_Alt_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_OTH_Ref_freq}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gnomAD_exome_OTH_Alt_freq}}{double COLUMN_DESCRIPTION}
#'}
#' @name pop_result
NULL

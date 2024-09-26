#' Find gene information from probe set
#'
#' The function `probes2genes()` is used to get gene name, symbol, entrezid id etc from probes id.
#'
#' @param affy_ids Affymetric probe IDs
#' @param organism Sample organism used for experiment ("human" or "rat")
#'
#' @return
#' A table of Ensembl gene symbol and other information
#' @export
#'
#' @examples
#' \dontrun{
#' affy_keys <- AnnotationDbi::keys(rat2302.db::rat2302.db, keytype = "PROBEID")
#' sample_id <- sample(affy_keys, 10)
#' probes2genes(affy_ids = sample_id, organism = "rat")
#' }
probes2genes <- function(affy_ids, organism = "rat") {
  if (identical(organism, "rat")) {
    gene_tab <- suppressMessages(AnnotationDbi::select(
      x = rat2302.db::rat2302.db,
      keys = affy_ids ,
      columns = c("SYMBOL","ENTREZID", "ENSEMBL","GENENAME")
    ))
  } else if(identical(organism, "human")) {
    gene_tab <- suppressMessages(AnnotationDbi::select(
      x = hgu133plus2.db::hgu133plus2.db,
      keys = affy_ids ,
      columns = c("SYMBOL","ENTREZID", "ENSEMBL","GENENAME")
    ))
  }
  gene_tab <- gene_tab[!duplicated(gene_tab$PROBEID), ]
  return(gene_tab)
}

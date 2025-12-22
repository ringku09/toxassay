#' Map Affymetrix probe IDs to gene annotations
#'
#' `probes2genes()` retrieves gene-level annotation for Affymetrix probe set IDs,
#' including gene symbols, Entrez IDs, Ensembl IDs, and gene names.
#'
#' @param affy_ids A character vector of Affymetrix probe set IDs.
#' @param organism A character string specifying the organism used in the experiment.
#'   Supported values are `"rat"` and `"human"`.
#'
#' @return
#' A data frame containing probe IDs and their corresponding gene annotations,
#' including gene symbol, Entrez ID, Ensembl ID, and gene name.
#'
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
  names(gene_tab) <- c("probe_id", "gene_symbol", "entrez_id", "ensembl_id", "gene_name")

    tolower(names(gene_tab))
  return(gene_tab)
}

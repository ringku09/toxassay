#' List of Compounds Available in the TG-GATES Database
#'
#' @description
#' The `tggates_compounds()` function retrieves a list of available compounds from
#' the [Open TG-GATEs](https://toxico.nibiohn.go.jp/english/) database. These compounds
#' can be used to download and process Affymetrix CEL images of samples.
#'
#' @section Open TG-GATEs Database:
#' The [TG-GATEs](https://toxico.nibiohn.go.jp/english/) database includes tests on 170
#' distinct compounds, primarily medicinal drugs, across both *in-vivo* and *in-vitro*
#' experimental setups. For *in-vivo* trials, male Sprague Dawley rats were used, and
#' the studies were organized into single-dose and repeated-dose protocols. Additionally,
#' the database includes two *in-vitro* studies utilizing primary hepatocytes from male
#' Sprague Dawley rats and human donors.
#'
#' @param species A character vector specifying the species, either `"Rat"` or `"Human"`.
#' @param data_type A character vector specifying the type of data, either `"in_vivo"` or `"in_vitro"`.
#' @param tissue A character vector specifying the organ, either `"Liver"` or `"Kidney"`.
#' @param dose_type A character vector specifying the dose type for *in-vivo* data, either `"Single"` or `"Repeat"`.
#'
#' @family Open TG-GATEs Data Download Helpers
#'
#' @return A vector of compound names available for the specified criteria.
#'
#' @export
#'
#' @examples
#' # Data available for compounds in Rat in-vivo liver Single data
#' tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#' # Data available for compounds in Rat in-vivo liver Repeat data
#' tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Repeat")
#' # Data available for compounds in Rat in-vitro liver Single data
#' tggates_compounds(species = "Rat", data_type = "in_vitro", tissue = "Liver")
#' # Data available for compounds in Rat in-vitro liver Repeat data
#' tggates_compounds(species = "Rat", data_type = "in_vitro", tissue = "Liver", dose_type = "Repeat")
#' # Data available for compounds in human in-vitro liver Single data
#' tggates_compounds(species = "Human", data_type = "in_vitro", tissue = "Liver")
#' # Data available for compounds in Rat in-vivo kidney single data
#' tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "kidney", dose_type = "single")
#' # Data available for compounds in Rat in-vivo kidney Repeat data
#' tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "kidney", dose_type = "Repeat")
tggates_compounds <- function(species = c("Rat", "Human"),
                              data_type = c("in_vivo", "in_vitro"),
                              tissue = c("Liver", "Kidney"),
                              dose_type = c("Single", "Repeat"),
                              error_call = caller_env()) {
  species <- block_fst(species)
  data_type <- gsub("\\W", "_", data_type)
  tissue <- block_fst(tissue)
  dose_type <- block_fst(dose_type)
  test_input(species, auto_input = TRUE)
  test_input(data_type, auto_input = TRUE)
  test_input(tissue, auto_input = TRUE)
  test_input(dose_type, auto_input = TRUE)
  tggates_data <- compounds_tggates %>%
    mutate(across(- c("compound_name", "compound_abbr"), ~ ifelse(is.na(.), 0, 1))) %>%
    mutate(across(where(is.numeric), ~ . == 1))
  if (identical(data_type, "in_vivo")) {
    query_column <- glue::glue("{species}_{data_type}_{tissue}_{dose_type}")
  } else {
    query_column <- glue::glue("{species}_{data_type}_{tissue}")
  }
  if (! query_column %in% names(tggates_data)[-(1:2)]) {
    cli_abort(c("Invalid data type parameter selection",
                "x" = "Data type {style_bold(col_red(backtick(query_column)))} is \\
                not available in the open TG-GATEs database.",
                "i" = " Please use one of the following valid parameter combinations: \\
                {style_italic(col_blue(backtick(names(tggates_data)[-(1:2)])))}.")
              , call = error_call)
  }
  tgp_com <- tggates_data %>%
    dplyr::filter(!!rlang::sym(query_column)) %>%
    dplyr::select(compound_name) %>%
    dplyr::pull() %>%
    sort()
  #cli_alert(glue("Data available for the data type {style_bold(col_green(backtick(query_column)))} of compounds:"))
  return(tgp_com)
}

#' List of Compounds Available in DrugMatrix Database
#'
#' @description
#' The `drugmatrix_compounds()` function retrieves a list of available compounds from
#' the [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) database. These compound names
#' can be used to download and process Affymetrix CEL images and metadata from
#' [GEO](https://www.ncbi.nlm.nih.gov/geo/).
#'
#' @section Available Data for Different Tissues:
#' The [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) database includes tests on
#' over 600 distinct compounds, conducted in both `in vivo` and `in vitro`
#' experimental setups. Male Sprague Dawley rats were used in both `in vivo` and `in vitro`
#' experiments, with `in vitro` studies specifically employing rat primary hepatocytes.
#' The available data for different tissues are:
#'    * Liver
#'    * Kidney
#'    * Heart
#'    * Hepatocytes (Rat Primary Hepatocytes)
#'
#' @param tissue Optional argument specifying the tissue type. Default is missing. See the section
#'   **Available Data for Different Tissues** for details.
#'
#'  One of:
#'    - `"Liver"`, for data from liver tissue.
#'    - `"Kidney"`, for data from kidney tissue.
#'    - `"Heart"`, for data from heart tissue.
#'    - `"Hepatocytes"`, for data from rat primary hepatocytes.
#'
#' @family DrugMatrix data download helpers
#'
#' @return
#' A data frame of boolean values if the parameter `tissue` is missing. The first column
#' contains available compound names. The subsequent columns represent different types of
#' available data. In the data frame:
#'  * `TRUE` indicates that data is available for the corresponding compound.
#'  * `FALSE` signifies that it is not.
#' For a given `tissue`, the output will be a vector containing the available compound names.
#'
#' @export
#'
#' @examples
#' # Data available for compounds in `Liver`
#' drugmatrix_compounds(tissue = "Liver")
#' # Data available for compounds in `Kidney`
#' drugmatrix_compounds(tissue = "Kidney")
#' # Data available for compounds in `Heart`
#' drugmatrix_compounds(tissue = "Heart")
#' # Data available for compounds in `Hepatocytes`
#' drugmatrix_compounds(tissue = "Hepatocytes")
drugmatrix_compounds <- function(tissue = c("Liver", "Kidney", "Heart", "Hepatocytes")) {
  dm_comps <- dm_metadata %>%
    dplyr::filter(Compound != "Control") %>%
    dplyr::select(c(Compound, Tissue)) %>%
    dplyr::distinct()
  dm_comps$values <- 1
  dm_data <- dm_comps %>%
    pivot_wider(names_from = Tissue, values_from = values, values_fill = 0) %>%
    rename("Hepatocytes (in vitro)" = Hepatocytes) %>%
    mutate(across(where(is.numeric), ~ . == 1))
    tissue <- block_fst(tissue)
    test_input(tissue, auto_input = TRUE)
    if (identical(tissue, "Liver")) {
      dm_com <- dm_data %>%
        dplyr::filter(Liver) %>%
        dplyr::select(Compound) %>%
        dplyr::pull( ) %>%
        sort()
    } else if (identical(tissue, "Kidney")) {
      dm_com <- dm_data %>%
        dplyr::filter(Kidney) %>%
        dplyr::select(Compound) %>%
        dplyr::pull() %>%
        sort()
    } else if (identical(tissue, "Heart")) {
      dm_com <- dm_data %>%
        dplyr::filter(Heart) %>%
        dplyr::select(Compound) %>%
        dplyr::pull() %>%
        sort()
    } else if (identical(tissue, "Hepatocytes")) {
      dm_com <- dm_data %>%
        dplyr::filter(`Hepatocytes (in vitro)`) %>%
        dplyr::select(Compound) %>%
        dplyr::pull() %>%
        sort()
    }
    #cli_alert(glue("Data available for the {style_bold(col_green(backtick(tissue)))} tissue of compounds:"))
  return(dm_com)
}

#' Collect Metadata for Compounds in the DrugMatrix Database
#'
#' @description
#' The `dm_metadata()` function retrieves metadata for query compounds from the
#' [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) database. It uses the compound name and tissue type to
#' obtain Affymetrix CEL images and their attribute data from the
#' [GEO](https://www.ncbi.nlm.nih.gov/geo/) database. Each CEL file represents raw expression
#' data of `probes` measured at various dose and time points, or combinations thereof, from chemical compound exposure.
#'
#' @inheritSection drugmatrix_compounds Available Data for Different Tissues
#'
#' @param compounds A vector of compound names to query.
#' @param tissue The tissue type used to generate perturbation data. Default is `Liver`.
#'   Refer to the **Available data for different tissue** section for more details.
#'
#' @return
#' A data frame containing the metadata.
#'
#' @export
#'
#' @examples
#' # Perturbation data of Ethanol from liver tissue
#' drugmatrix_dictionary("Ethanol", "Liver")
#' # Perturbation data of Ethanol from Kidney tissue
#' drugmatrix_dictionary("Ethanol", "Kidney")
drugmatrix_dictionary <- function(compounds, tissue = c("Liver", "Kidney", "Heart", "Hepatocytes")) {
  tissue <- block_fst(tissue)
  test_input(tissue, auto_input = TRUE)
  vehicle <- unique(dm_metadata$Vehicle[dm_metadata$Compound %in% compounds])

  pert_meta <- dm_metadata[dm_metadata$Compound %in% compounds & dm_metadata$Tissue == tissue, ]
  ctrl_meta <- dm_metadata[dm_metadata$Compound == "Control"
                           & dm_metadata$Tissue == tissue
                           & dm_metadata$Vehicle %in% vehicle
                           & dm_metadata$Time %in% unique(pert_meta$Time), ]
  attr_data <- dplyr::bind_rows(pert_meta, ctrl_meta)
  return(attr_data)
}

#' Download Link for Data in the Open TG-GATEs Database
#'
#' @description
#' The `tggates_link()` function retrieves the download link for perturbation data of compounds along with metadata from the [TG-GATEs](https://toxico.nibiohn.go.jp/english/) database.
#'
#' @param compound An available compound name.
#' @inheritSection tggates_compounds Open TG-GATEs Database
#'
#' @family Open TG-GATEs Data Download Helpers
#'
#' @return A download link.
#' @export
#'
#' @examples
#' # Download link for downloading perturbation data of Acetaminophen for different experimental settings
#' tggates_link(compound = "acetaminophen", species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "single")
#' tggates_link(compound = "acetaminophen", species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Repeat")
#' tggates_link(compound = "acetaminophen", species = "Rat", data_type = "in_vitro", tissue = "Liver", dose_type = "Single")
#' tggates_link(compound = "acetaminophen", species = "Human", data_type = "in_vitro", tissue = "Liver", dose_type = "Single")
tggates_link <- function(compound,
                         species = c("Rat", "Human"),
                         data_type = c("in_vivo", "in_vitro"),
                         tissue = c("Liver", "Kidney"),
                         dose_type = c("Single", "Repeat")) {
  species <- block_fst(species)
  data_type <- gsub("\\W", "_", data_type)
  tissue <- block_fst(tissue)
  dose_type <- block_fst(dose_type)
  test_input(species, auto_input = TRUE)
  test_input(data_type, auto_input = TRUE)
  test_input(tissue, auto_input = TRUE)
  test_input(dose_type, auto_input = TRUE)
  is_compound(compounds = compound,
              database = "tggates",
              tissue = tissue,
              species = species,
              data_type = data_type,
              dose_type = dose_type)
  com_info <-  glue::glue("{compound}.{species}.{data_type}.{tissue}.{dose_type}.zip")
  tggates_url <- "https://dbarchive.biosciencedbc.jp/data/open-tggates/LATEST"
  com_link <- glue::glue("{tggates_url}/{species}/{data_type}/{tissue}/{dose_type}/{com_info}")
  #com_link <- as.character(com_link)
  return(com_link)
}

#' Download Links for Data in the DrugMatrix Database
#'
#' @description
#' The `drugmatrix_link()` function retrieves download links for the perturbation data
#' of specified compounds from the [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) database.
#'
#' @inheritSection drugmatrix_compounds  Available Data for Different Tissues
#'
#' @inheritParams drugmatrix_dictionary
#'
#' @family DrugMatrix data download helpers
#'
#' @return
#' A vector of download links for the specified compounds and tissue.
#'
#' @export
#'
#' @examples
#' # Downloading link for perturbation data of Ethanol from Liver tissue
#' drugmatrix_link(compounds= "Ethanol", tissue = "Liver")
#' # Downloading link for perturbation data of Ethanol and Esrriol from Liver tissue
#' drugmatrix_link(c("Ethanol", "Estriol"))
drugmatrix_link <- function(compounds, tissue = c("Liver", "Kidney", "Heart", "Hepatocytes")) {
  tissue <- block_fst(tissue)
  test_input(tissue, auto_input = TRUE)
  is_compound(compounds, database = "drugmatrix", tissue = tissue)
  meta <- drugmatrix_dictionary(compounds = compounds, tissue = tissue)
  cel_links <- mapply(
    function(gsm_id, sample_id)
      glue::glue("https://www.ncbi.nlm.nih.gov/geo/download/?acc={gsm_id}&format=file&file={gsm_id}%5F{sample_id}%2ECEL%2Egz"),
    gsm_id = meta$Accession,
    meta$Sample_id
  )
  return(cel_links)
}

#' Download Link for Data in the CTD Database
#'
#' @description
#' The `ctd_link()` function retrieves the download link for perturbation data of compounds
#' from the [CTD (Comparative Toxicogenomics Database)](http://ctdbase.org/).
#'
#' @param compounds A character vector specifying the names of the compounds.
#' @param association A character string specifying the type of association to retrieve (e.g., "genes_curated"). Default is `"cgixns"`.
#' @param action_type A character string specifying the type of action (e.g., "ANY"). Default is `"ANY"`.
#'
#' @family CTD data download helpers
#'
#' @return A character vector containing the download links for the specified compounds.
#' @export
#'
#' @examples
#' # Downloading link for association data of Ethanol and curated genes
#' ctd_link(compounds = "Ethanol", association = "genes_curated")
#' # Downloading link for association data of all type between Ethanol and Esrriol, and genes
#' ctd_link(compounds = c("Ethanol", "Estriol"), association = "cgixns")
ctd_link <- function(compounds,
                     association = c("cgixns", "genes_curated", "genes_inferred",
                                     "diseases", "diseases_curated", "diseases_inferred"),
                     action_type = "ANY") {
  compounds <- gsub(" ", "%20", compounds)
  ctd_compounds <- paste(compounds, collapse = "|")
  test_input(association, auto_input = TRUE)
  if (identical(association, "cgixns")) {
    if (is.null(action_type)) {
      cli_abort(c(
        ".arg(action_type) must provided for association of type `cgixns`",
        "x" = "There is no data for {style_bold(col_red(backtick(compounds[wrong_com])))} compound{?s} in {style_bold(col_br_red(backtick(tissue)))} tissue.",
        "i" = "To get available compounds in {style_bold(col_br_red(backtick(tissue)))} tissue, call the function {style_italic(col_blue(backtick(avail_com)))}."
      ))
    }
    query_link <- glue::glue("https://ctdbase.org/tools/batchQuery.go?inputType=chem&inputTerms={ctd_compounds}&report={association}&actionTypes={action_type}&format=tsv")
  } else {
    query_link <- glue::glue("https://ctdbase.org/tools/batchQuery.go?inputType=chem&inputTerms={ctd_compounds}&report={association}&format=tsv")
  }
  return(query_link)
}

#' Download Microarray and Metadata for Query Compounds from the TG-GATEs Database
#'
#' @description
#' The `tggates_rawdata()` function downloads Affymetrix CEL files and associated attribute data,
#' including a data dictionary, clinical, and biochemical information, for a specific compound
#' from the [TG-GATEs](https://toxico.nibiohn.go.jp/english/) database. The data will be
#' organized in a directory named **compound_name.data_type** in the specified output location.
#' Within this directory, a sub-folder named **celfiles** will store all CEL files, and a .tsv
#' file containing attribute information will also be created.
#'
#' The function ensures that gene expression data, measured at various doses and time points,
#' are systematically stored. The download speed will depend on your internet connection.
#'
#' @inheritParams tggates_link
#' @param output_dir Output directory for downloading data. If not specified,
#'  the temporary path of the current R session will be used.
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A folder containing `CEL` files and a `tsv` metadata file.
#' @export
#'
#' @examples
#' \dontrun{
#' tggates_rawdata(compound = "2-nitrofluorene", species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#' tggates_rawdata(compound = "N-nitrosomorpholine", species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#' }
tggates_rawdata <- function(compound,
                            species = c("Rat", "Human"),
                            data_type = c("in_vivo", "in_vitro"),
                            tissue = c("Liver", "Kidney"),
                            dose_type = c("Single", "Repeat"),
                            output_dir = missing_arg()) {
  # change the download need for base `download.file` function timeout since the files are big
  # opts <- options()
  # options(timeout=1000)
  # on.exit(options(opts))
  # options(timeout = max(1000, getOption("timeout")))
  output_dir <- destination(output_dir)
  data_link <- tggates_link(compound = compound,
                                  species = species,
                                  data_type = data_type,
                                  tissue = tissue,
                                  dose_type = dose_type)
  `Open_TG-GATEs` <- "https://dbarchive.biosciencedbc.jp/data/open-tggates/"
  check_internet(`Open_TG-GATEs`)
  temp_comp <- gsub("^.*/", "", data_link)
  sp_comp <- unlist(strsplit(temp_comp, "\\."))
  dt_type <- glue_collapse(sp_comp[-c(1, length(sp_comp))],"-")
  out_path <- glue::glue("{output_dir}/{temp_comp}")
  start_time <- Sys.time()
  #utils::download.file(comp_link, out_path, quiet = TRUE)
  #curl::curl_download(comp_link, out_path)
  cli_alert_info("Downloading open TG-GATEs perturbation data for {style_bold(col_red(backtick(sp_comp[1])))} ({col_br_red(dt_type)}):")
  httr::GET(data_link,httr::write_disk(out_path, overwrite = TRUE), httr::progress())
  end_time <- Sys.time()
  req_time <- difftime(end_time, start_time, units = "secs")[[1]]
  file_size <- file.info(out_path)$size
  zip_file <- list.files(path = output_dir, pattern = temp_comp, full.names = TRUE)
  utils::unzip(zip_file, exdir = output_dir)
  rm_zip <- file.remove(zip_file)
  #sp_comp <- unlist(strsplit(temp_comp, "\\."))
  #dt_type <- glue_collapse(sp_comp[-c(1, length(sp_comp))],"-")
  if(rm_zip) {
    cli_alert_success(c("Downloaded {prettyunits::pretty_bytes(file_size)}",
                             " in {prettyunits::pretty_sec(as.numeric(req_time))}"))
  }else {
    cli_alert_success("{style_bold(col_red(backtick(sp_comp[1])))} ({col_br_red(dt_type)}) data download has unsuccessful")
  }
  path <- gsub(".zip$","" , out_path)
  attr_file <- paste(path, "Attribute.tsv", sep = "/")
  utils::read.table(file = attr_file, sep = '\t', header = TRUE) %>%
    dplyr::mutate(database = "tggates") %>%
    readr::write_tsv(attr_file)
  return(path)
}

#' Download Microarray and Metadata for Query Compounds from the DrugMatrix Database
#'
#' @description
#' The `drugmatrix_rawdata()` function downloads raw Affymetrix CEL files and their associated
#' attribute files (data dictionary) for a specified compound from the
#' [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) database. A folder named
#' **celfiles** will be created in your specified output directory to store the CEL image
#' files, and an accompanying TSV file containing metadata will also be generated in the
#' parent folder.
#'
#' This function retrieves Affymetrix CEL image files for specific compounds from the
#' DrugMatrix database, storing all downloaded data in a user-selected directory. Since gene
#' expression is measured at various doses and time points or combinations thereof, a
#' collection of CEL files will be downloaded into a **celfiles** folder. Additionally, a
#' TSV file containing comprehensive metadata will be saved in the parent directory. Note
#' that download speed depends on your internet connection.
#'
#' @inheritParams drugmatrix_dictionary
#' @inheritParams tggates_rawdata
#'
#' @family DrugMatrix data download helpers
#'
#' @return A folder containing `CEL` files and a `TSV` metadata file.
#' @export
#'
#' @examples
#' compounds <- "Ethanol"
#' \dontrun{
#' drugmatrix_rawdata(compounds, "Liver")
#' }
drugmatrix_rawdata <- function(compounds,
                       tissue = c("liver", "kidney", "heart", "hepatocytes"),
                        output_dir = missing_arg()) {
  output_dir <- destination(output_dir)
  output_path <- glue::glue("{output_dir}/celfiles")
  dir.create(output_path)
  comp_link <- drugmatrix_link(compounds = compounds, tissue = tissue)
  DrugMatrix <- "https://www.ncbi.nlm.nih.gov/geo/download/"
  check_internet(DrugMatrix)
  cli_alert_info("Downloading DrugMatrix perturbation data for {style_bold(col_red(backtick(compounds)))} in {col_green(tissue)} tissue:")
  start_time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = length(comp_link), style = 3)
  for (i in seq_along(comp_link)) {
    out_path <- glue::glue("{output_dir}/celfiles/{names(comp_link[i])}.CEL")
    # utils::download.file(url = comp_link[i], destfile = out_path, quiet = TRUE)
    httr::GET(comp_link[i],httr::write_disk(out_path, overwrite = TRUE))
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  end_time <- Sys.time()
  req_time <- difftime(end_time, start_time, units = "secs")[[1]]
  files <- list.files(path = output_dir, full.names = TRUE, recursive = TRUE)
  file_info <- file.info(files)
  total_size <- sum(file_info$size, na.rm = TRUE)
  if(i == length(comp_link)) {
    cli_alert_success(c("Downloaded {prettyunits::pretty_bytes(total_size)} ",
                        " in {prettyunits::pretty_sec(as.numeric(req_time))}"), wrap = TRUE)
  }else {
    cli_alert_success("Data download has unsuccessful. Please try again or check your internet connection")
  }
  attr_data <- drugmatrix_dictionary(compounds, tissue) %>%
    dplyr::mutate(database = "drugmatrix")
  readr::write_tsv(attr_data, glue::glue("{output_dir}/Attribute.tsv"))
  return(output_dir)
}

#' Convert Affymetrix CEL Files to Gene Expression Matrix
#'
#' @description
#' The `cel2exprs()` function extracts gene expression data from CEL files, outputting
#' an expression matrix. The dimensions of the matrix depend on the number of CEL files
#' and the probes specified in the Chip Definition File (CDF).
#'
#' This function utilizes core functions [`mas5`][affy::mas5()] and [`rma`][affy::rma()]
#' from the affy package to quantify mRNA intensity of probes/genes. The expression
#' measurements are provided on a log base 2 scale.
#'
#' @param cel_files A vector containing the full paths of CEL files.
#' @param organism The sample organism used in the experiment, either `"human"` or `"rat"`.
#' @param method The normalization method to use, either `"rma"` or `"mas5"`. Default is `"rma"`.
#' @param ... Additional options for processing CEL files, passed on to
#'   [rma][affy::rma()] or [`mas5`][affy::mas5()], such as `normalize`, `background`, etc.
#'
#' @family Gene Expression Data Processing
#'
#' @seealso
#'    [affy::mas5()] for Affymetrix version 5 (MAS5) model,
#'    [affy::rma()] for Robust Multi-array Average (RMA) model
#'
#' @return
#' A matrix of raw gene expression data:
#'   * Rows represent probes
#'   * Columns represent sample CEL IDs
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- tggates_rawdata(compound = "2-nitrofluorene")
#' cel_file <- file.path(file_path, "celfiles")
#' cel_files <- list.files(cel_file, full.names = TRUE)
#' expr_data <- cel2exprs(cel_files)
#' head(expr_data)
#' }
cel2exprs <- function(cel_files, organism = "rat", method = "rma", ...)
  {
  if (identical(organism, "human")) {
    CDF = "hgu133plus2cdf"
  } else if(identical(organism, "rat")) {
    CDF = "rat2302cdf"
  } else {
    cli_abort(c(
      "Wrong organism name",
      "x" = "The name of {.var organism} ({style_bold(col_blue(backtick(organism)))}) you provided is not recognized",
      "i" = "Please chosse either `human` or `rat` instead"
      ))
  }
  files <- tools::file_ext(cel_files)
  cel_path <- cel_files[which(files == "CEL")]
  if (!all(files == "CEL")) {
    cli_warn(c(
      "Some of the paths you provided do not contain any CEL files.",
      i = "The following paths have been ignored:",
      format_bullets_raw(rlang::set_names(cel_files[!files == "CEL"], "x"))
      ))
  }
  read_cel <- affy::ReadAffy(filenames = cel_path, cdfname = CDF)
  if (identical(method, "mas5")) {
    expr_meas <- affy::mas5(read_cel, verbose = FALSE, ...)
    expr <- Biobase::exprs(expr_meas)
    expr <- log2(expr)
  } else if (identical(method, "rma")) {
    expr_meas <- affy::rma(read_cel, verbose = FALSE, ...)
    expr <- Biobase::exprs(expr_meas)
  }
  return(expr)
}

#' Processing Gene Expression and Metadata of a Compound in the TG-GATEs Database
#'
#' @description
#' The `process_celfile()` function generates gene expression data and associated metadata for samples.
#' This data is derived from CEL files and attribute files, respectively.
#'
#' @param compound_path Full path to the compound directory, which should include
#'  a folder named `celfile` containing CEL files, as well as a `tsv` attribute file.
#' @param fc A boolean indicating the format of the output gene expression data.
#'   * `TRUE` (default): The output will be in fold-change (FC) format.
#'   * `FALSE`: The output will be raw gene expression data.
#' @param multicore Used for parallel computing.
#'   * `FALSE` (default): No parallelization will be used.
#'   * `TRUE`: Computation will employ parallel processing, leaving one core for system use.
#'   * Integer: Specifies the number of workers for parallel processing.
#' @inheritParams tggates_rawdata
#' @inheritParams cel2exprs
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A list containing:
#'   * `gexpr_data`: A matrix of gene expression values, with probes in rows and samples (CEL IDs) in columns.
#'   * `meta_data`: A data frame containing metadata about the samples and other related information.
#' @export
#'
#' @examples
#' compound <- "2NF"
#' \dontrun{
#' file_path <- tggates_rawdata(compound = "2-nitrofluorene")
#' proc_expr <- process_celfile(file_path)
#' head(proc_expr$gexpr_data)
#' head(proc_expr$meta_data)
#' }
process_celfile <- function(compound_path,
                            fc = TRUE,
                            method = "rma",
                            multicore = FALSE,
                            output_dir = missing_arg(),
                            ...) {
  output_dir <- destination(output_dir)
  start_time <- Sys.time()
  sp_path <- unlist(strsplit(compound_path, "\\."))
  sp_comp <- unlist(strsplit(sp_path[1], "/"))
  comp_nm <- sp_comp[length(sp_comp)]
  comp_typ <- glue_collapse(sp_path[-1],"-")
  start_msg <- "Processing {style_bold(col_red(backtick(comp_nm)))} ({col_br_red(comp_typ)}) data........."
  cat(glue(start_msg))
  attr_file <- paste(compound_path, "Attribute.tsv", sep = "/")
  #comp_dict <- readr::read_tsv(file = attr_file, show_col_types = FALSE)  # not work for rat data
  comp_dict <- utils::read.table(file = attr_file, sep = '\t', header = TRUE) # not work for human data
  comp_dict <- comp_dict %>%
    dplyr::rename_with(~ tolower(gsub(" ", "_", .x, fixed = TRUE))) %>%
    dplyr::rename(
      compound_abbr = "compound.abbr.",
      time_level = "sacri_period",
      sample_id = "individual_id"
    ) %>%
    dplyr::filter(.data$barcode != "No ChipData") %>%
    dplyr::mutate_at(c("compound_name", "compound_abbr", "dose_level", "time_level"), as.factor) %>%
    dplyr::relocate(
      c("compound_name", "compound_abbr", "dose_level", "time_level",
        "sample_id", "arr_design", "species", "test_type", "organ_id", "sin_rep_type"),
      .after = "barcode")
  chip <- unique(comp_dict$arr_design)
  if(identical(chip, "Rat230_2")) {
    organism = "rat"
  }else if(identical(chip, "HG-U133_Plus_2")) {
    organism = "human"
  }
  cel_file <- paste(compound_path, "celfiles", sep = "/")
  file_path <- list.files(path = cel_file, pattern = "*.CEL", full.names = TRUE)
  if(is.logical(multicore)) {
    if(multicore ) {
      multicore <- start_parallel(multicore)
      stop_cluster <- TRUE
    } else {
      multicore <- stop_cluster <- FALSE
    }
  }else {
    stop_cluster <- if(inherits(multicore, "cluster")) FALSE else TRUE
    multicore <- start_parallel(multicore)
  }
  on.exit(if(multicore & stop_cluster)
    stop_parallel(attr(multicore, "cluster")))
  `%DO%` <- if(multicore) foreach::`%dopar%` else foreach::`%do%`
  if(!multicore) {
    # cli_alert_warning(c("Processing without parallel computing will take much time. ",
    #                   "Please use `multicore = TRUE` to parallelize your task."))
    expr_data <- cel2exprs(file_path, organism = organism, method = method, ...)
  }else {
    ii <- attr(multicore, "cores")
    # Set condition if ii>length(file_path), then use length(file_path) cores only.
    par_idx <- parallel::splitIndices(length(file_path), ii)
    fun_call <- as.call(c(list(quote(foreach::foreach), i = seq_len(ii)), .combine ="cbind"))
    .fun <- eval(fun_call)
    #cli_alert_warning(c("Processing will takes time."))
    expr_data <- .fun %DO% cel2exprs(file_path[par_idx[[i]]], organism = organism, method = method, ...)
  }
  colnames(expr_data) <- gsub(".CEL", "", colnames(expr_data))
  timepnt <- levels(comp_dict$time_level)
  temp_dict <- comp_dict %>%
    dplyr::filter(.data$dose_level != "Control")
  if(fc) {
    temp_expr <- list()
    for(i in 1:length(timepnt)) {
      const_dict <- comp_dict %>%
        dplyr::filter(.data$dose_level == "Control" & .data$time_level == timepnt[i])
      treat_dict <- comp_dict %>%
        dplyr::filter(.data$dose_level != "Control" & .data$time_level == timepnt[i])
      const_expr <- expr_data[, colnames(expr_data) %in% const_dict$barcode]
      treat_expr <- expr_data[, colnames(expr_data) %in% treat_dict$barcode]
      const_med <- apply(const_expr, 1, stats::median)
      temp_expr[[i]] <- treat_expr - const_med
    }
    raw_expr <- do.call(cbind, temp_expr)
    match_col <-  match(temp_dict$barcode , colnames(raw_expr))
    match_col <- match_col[!is.na(match_col)]
    exprs_data <-raw_expr[, match_col]
    temp_dict <- temp_dict %>%
      mutate(fc = TRUE, .before = "arr_design")
  }else {
    match_col <-  match(temp_dict$barcode, colnames(expr_data))
    match_col <- match_col[!is.na(match_col)]
    exprs_data <- expr_data[, match_col]
  }
  out_data <- list(exprs_data, temp_dict)
  names(out_data) <- c("gexpr_data", "meta_data")
  end_time <- Sys.time()
  if(!multicore) {
    req_time <- difftime(end_time, start_time, units = "secs")[[1]]
    rm_msg(start_msg)
    cli_alert_success(c("Processing TG-GATEs perturbation data for {style_bold(col_red(backtick(comp_nm)))}",
                    "({col_br_red(comp_typ)}) has been completed in ",
                    "{prettyunits::pretty_sec(as.numeric(req_time))}."))
  }else {
    req_time <- difftime(end_time, start_time, units = "secs")[[1]]
    rm_msg(start_msg)
    cli_alert_success(c("Processing TG-GATEs perturbation data for {style_bold(col_red(backtick(comp_nm)))}",
                        "({col_br_red(comp_typ)}) has been completed ",
                        "in {prettyunits::pretty_sec(as.numeric(req_time))} using {ii} worker{?s}."))
  }
  return(out_data)
}

#' Process Compound Data from Open TG-GATEs
#'
#' @description
#' The `process_tggates()` function processes a set of compounds from the Open TG-GATEs
#' database. It handles all CEL and attribute files for each compound and generates two
#' CSV files: `expression.csv` and `metadata.csv`. If the `store` flag is set to `TRUE`,
#' these CSV files will be saved in the specified output directory.
#'
#' @param files A vector of paths to compound directories, each containing a `celfile` folder
#'   and a `.tsv` attribute file.
#' @inheritParams process_celfile
#' @param store Logical. If `TRUE` (the default), the resulting data will be saved as
#'   two CSV files (gene expression and metadata) in the output directory.
#' @param ... Additional arguments passed to `process_celfile`.
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A list containing:
#'   * `expression`: A matrix of gene expression values, with probes in rows and samples (CEL IDs) in columns.
#'   * `metadata`: A data frame with attributes about samples and clinical information of compounds.
#' @export
#'
#' @examples
#' \dontrun{
#' compound_list <- list("2-nitrofluorene", "N-nitrosomorpholine")
#' files <- lapply(compound_list, tggates_rawdata, species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#' tgx_data <- process_tggates(files, store = TRUE)
#' head(tgx_data$expression)
#' head(tgx_data$metadata)
#' }
process_tggates <- function(files,
                            fc = TRUE,
                            method = "rma",
                            multicore = FALSE,
                            store = FALSE,
                            output_dir = missing_arg(),
                            ...) {
  files <- as.list(files)
  exprs_list <- lapply(files, process_celfile, fc = fc,
                       multicore = multicore, method = method, output_dir = output_dir, ...)
  expr_mat <- do.call(cbind, lapply(exprs_list, "[[",1))
  dict_df <- do.call(dplyr::bind_rows, lapply(exprs_list, "[[",2))
  if (store) {
    output_dir <- destination(output_dir)
    expr_mat %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    dict_df %>%
      readr::write_csv(file = paste(output_dir, "metadata.csv", sep = "/"))
  }
  tgx_data <- list(expr_mat, dict_df)
  names(tgx_data) <- c("expression", "metadata")
  return(tgx_data)
}

#' Process Gene Expression and Metadata from DrugMatrix Database
#'
#' @description
#' The `process_drugmatrix()` function processes gene expression data and metadata for samples
#' from CEL files and attribute files in the DrugMatrix database.
#'
#' @param rawdata_path A string specifying the full path to the raw data directory of compounds, which should contain a folder named `celfile` and a `tsv` attribute file.
#' @inheritParams process_tggates
#'
#' @family Open TG-GATEs data download helpers
#'
#' @return A list with the following components:
#'   * `expression`: A matrix of gene expression values with probes in rows and samples (CEL IDs) in columns.
#'   * `metadata`: A data frame containing metadata about the samples.
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- drugmatrix_rawdata(compounds = "Ethanol", "Liver")
#' proc_dm <- process_drugmatrix(file_path)
#' head(proc_dm$expression)
#' head(proc_dm$metadata)
#' }
process_drugmatrix <- function(rawdata_path,
                       fc = TRUE,
                       method = "rma",
                       multicore = FALSE,
                       store = FALSE,
                       output_dir = missing_arg(),
                       ...) {
  comp_dict <- utils::read.table(file = glue("{rawdata_path}/Attribute.tsv"), sep = '\t', header = TRUE) # not work for human data
  comp_nm <- unique(comp_dict$Compound[comp_dict$Compound != "Control"])
  comp_typ <- unique(comp_dict$Tissue)
  cli_alert_info("Processing DrugMatrix perturbation data for {style_bold(col_red(backtick(comp_nm)))} in {col_green(comp_typ)} tissue......")
  #comp_dict <- readr::read_tsv(file = attr_file, show_col_types = FALSE)  # not work for rat data
  cel_file <- paste(rawdata_path, "celfiles", sep = "/")
  file_path <- list.files(path = cel_file, pattern = "*.CEL", full.names = TRUE)
  if(is.logical(multicore)) {
    if(multicore) {
      multicore <- start_parallel(multicore)
      stop_cluster <- TRUE
    } else {
      multicore <- stop_cluster <- FALSE
    }
  }else {
    stop_cluster <- if(inherits(multicore, "cluster")) FALSE else TRUE
    multicore <- start_parallel(multicore)
  }
  on.exit(if(multicore & stop_cluster)
    stop_parallel(attr(multicore, "cluster")))
  `%DO%` <- if(multicore) foreach::`%dopar%` else foreach::`%do%`
  start_time <- Sys.time()
  if(!multicore) {
    # cli_alert_warning(c("Processing without parallel computing will take much time. ",
    #                   "Please use `multicore = TRUE` to parallelize your task."))
    expr_data <- cel2exprs(file_path, organism = "rat", method = method, ...)
  }else {
    ii <- attr(multicore, "cores")
    # Set condition if ii>length(file_path), then use length(file_path) cores only.
    par_idx <- parallel::splitIndices(length(file_path), ii)
    fun_call <- as.call(c(list(quote(foreach::foreach), i = seq_len(ii)), .combine ="cbind"))
    .fun <- eval(fun_call)
    #cli_alert_warning(c("Processing will takes time."))
    expr_data <- .fun %DO% cel2exprs(file_path[par_idx[[i]]], organism = "rat", method = method, ...)
  }
  colnames(expr_data) <- gsub(".CEL", "", colnames(expr_data))
  timepnt <- unique(comp_dict$Time)
  temp_dict <- comp_dict %>%
    dplyr::filter(.data$Compound != "Control")
  if(fc) {
    temp_expr <- list()
    for(i in 1:length(timepnt)) {
      const_dict <- comp_dict %>%
        dplyr::filter(.data$Compound == "Control" & .data$Time == timepnt[i])
      treat_dict <- comp_dict %>%
        dplyr::filter(.data$Compound != "Control" & .data$Time == timepnt[i])
      const_expr <- expr_data[, colnames(expr_data) %in% const_dict$Accession]
      treat_expr <- expr_data[, colnames(expr_data) %in% treat_dict$Accession]
      const_med <- apply(const_expr, 1, stats::median)
      temp_expr[[i]] <- treat_expr - const_med
    }
    raw_expr <- do.call(cbind, temp_expr)
    match_col <-  match(temp_dict$Accession , colnames(raw_expr))
    match_col <- match_col[!is.na(match_col)]
    exprs_data <-raw_expr[, match_col]
    temp_dict <- temp_dict %>%
      mutate(fc = TRUE)
  }else {
    match_col <-  match(temp_dict$barcode, colnames(expr_data))
    match_col <- match_col[!is.na(match_col)]
    exprs_data <- expr_data[, match_col]
  }
  end_time <- Sys.time()
  if(!multicore) {
    req_time <- difftime(end_time, start_time, units = "secs")[[1]]
    # rm_msg(start_msg)
    cli_alert_success(c("Data processing has been completed in ",
                        "{prettyunits::pretty_sec(as.numeric(req_time))}."))
  }else {
    req_time <- difftime(end_time, start_time, units = "secs")[[1]]
    # rm_msg(start_msg)
    cli_alert_success(c("Data processing has been completed ",
                        "in {prettyunits::pretty_sec(as.numeric(req_time))} using {ii} worker{?s}."))
  }
  if (store) {
    output_dir <- destination(output_dir)
    exprs_data %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    temp_dict %>%
      readr::write_csv(file = paste(output_dir, "metadata.csv", sep = "/"))
  }
  dm_dict <- temp_dict %>%
    dplyr::rename_with(tolower, .cols = everything()) %>%
    dplyr::rename(barcode = accession,
           compound_name = compound,
           time_level = time,
           organ_id = tissue) %>%
    dplyr::mutate(arr_design = "Rat230_2")
  dm_data <- list(exprs_data, dm_dict)
  names(dm_data) <- c("expression", "metadata")
  return(dm_data)
}

#' Retrieve and Process TG-GATES Data
#'
#' @description
#' The `get_tggates()` function retrieves raw data for specified compounds from the TG-GATES database and processes it into gene expression values and metadata.
#'
#' @param compounds A character vector of compound names to retrieve from the TG-GATES database.
#' @param data_type A character string specifying the type of data to retrieve. See `tggates_rawdata` for available options.
#' @param fc A logical value indicating whether to apply fold-change calculations. Default is `TRUE`.
#' @param method A character string specifying the preprocessing method to use. Default is `"rma"`.
#' @param multicore A logical value indicating whether to use multiple cores for processing. Default is `FALSE`.
#' @param store A logical value indicating whether to store the processed data. Default is `FALSE`.
#' @param output_dir A character string specifying the directory to store output files. Default is `missing_arg()`.
#' @param ... Additional arguments passed to the `process_tggates` function.
#'
#' @return A list containing processed data from the TG-GATES database.
#' @export
#'
#' @examples
#' compounds <- c("2-nitrofluorene", "N-nitrosomorpholine")
#' \dontrun{
#' data <- get_tggates(compounds, species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#' head(data$expression)
#' head(data$metadata)
#' }
get_tggates <- function(compounds,
                        species = c("Rat", "Human"),
                        data_type = c("in_vivo", "in_vitro"),
                        tissue = c("Liver", "Kidney"),
                        dose_type = c("Single", "Repeat"),
                        fc = TRUE,
                        method = "rma",
                        multicore = FALSE,
                        store = FALSE,
                        output_dir = missing_arg(),
                        ...) {
  compounds <- as.list(compounds)
  raw_files <- lapply(
    compounds,
    tggates_rawdata,
    species = species,
    data_type = data_type,
    tissue = tissue,
    dose_type = dose_type,
    output_dir = output_dir)
  tg_data <- process_tggates(
    raw_files,
    fc = fc,
    multicore= multicore,
    method = method,
    output_dir = output_dir,
    store = store,
    ...)
  return(tg_data)
}


#' Retrieve and Process DrugMatrix Data
#'
#' @description
#' The `get_dm()` function retrieves raw data files for specified compounds and tissue types
#' from the DrugMatrix database and processes them to produce gene expression data and metadata.
#'
#' @param compounds A character vector specifying the names of the compounds to retrieve.
#' @param tissue A character string specifying the tissue type. Default is `"Liver"`.
#' @param fc A logical value indicating whether to apply fold change. Default is `TRUE`.
#' @param method A character string specifying the method for processing the CEL files. Default is `"rma"`.
#' @param store A logical value indicating whether to store the processed data. Default is `FALSE`.
#' @param multicore A logical value indicating whether to use multiple cores for processing. Default is `FALSE`.
#' @param output_dir A character string specifying the directory to store output files. Default is `missing_arg()`.
#' @param ... Additional arguments passed to the `process_drugmatrix` function.
#'
#' @return A list with the following components:
#'   * `expression`: A matrix of gene expression values with probes in rows and samples (CEL IDs) in columns.
#'   * `metadata`: A data frame containing metadata about the samples.
#' @export
#'
#' @examples
#' \dontrun{
#' dm_data <- get_drugmatrix(compounds = "Ethanol", tissue = "Liver")
#' head(dm_data$expression)
#' head(dm_data$metadata)
#' }
get_drugmatrix <- function(compounds,
                   tissue = "Liver",
                   fc = TRUE,
                   method = "rma",
                   store = FALSE,
                   multicore = FALSE,
                   output_dir = missing_arg(),
                   ...) {
  raw_files <- drugmatrix_rawdata(compounds = compounds, tissue = tissue, output_dir = output_dir)
  dm_data <- process_drugmatrix(raw_files,
                        fc = fc,
                        method = method,
                        multicore = multicore,
                        store = store,
                        output_dir = output_dir,
                        ...)
  return(dm_data)
}

#' Retrieve CTD Data for Compounds
#'
#' @description
#' The `get_ctd()` function retrieves data for specified compounds from the
#' [Comparative Toxicogenomics Database (CTD)](http://ctdbase.org/), including gene and disease associations.
#'
#' @param compounds A character vector specifying the names of the compounds.
#' @param gene_parm A character string specifying the gene association parameter. Default is `"cgixns"`.
#' @param disease_parm A character string specifying the disease association parameter. Default is `"diseases"`.
#' @param action_type A character string specifying the type of action (e.g., "ANY"). Default is `"ANY"`.
#' @param output_dir A character string specifying the directory to store output files. Default is `missing_arg()`.
#'
#' @family CTD data retrieval helpers
#'
#' @return A list with the following components:
#'   * `compound_gene`: A data frame containing the gene associations for the specified compounds.
#'   * `compound_disease`: A data frame containing the disease associations for the specified compounds.
#' @export
#'
#' @examples
#' \dontrun{
#' ctd_data <- get_ctd(compounds = c("Ethanol", "Estriol"), gene_parm = "cgixns", disease_parm = "diseases")
#' head(ctd_data$compound_gene)
#' head(ctd_data$compound_disease)
#' }
get_ctd <- function(compounds,
                    gene_parm = "cgixns",
                    disease_parm = "diseases",
                    action_type = "ANY",
                    store = FALSE,
                    output_dir = missing_arg()) {
  CTD <- "https://ctdbase.org/"
  check_internet(CTD)
  com_gene <- ctd_link(compounds = compounds,
                  association = gene_parm,
                  action_type = action_type)
  com_dis <- ctd_link(compounds = compounds,
                  association = disease_parm,
                  action_type = action_type)
  cli_alert_info("Downloading CTD data for `compound-gene` and `Compound-disease`
                 association using parameter {style_bold(col_red(backtick(gene_parm)))} and
                 {style_bold(col_red(backtick(disease_parm)))}, respectively", wrap = TRUE)
  # gene_data <- httr::GET(com_gene,
  #                        httr::write_disk(gene_path, overwrite = TRUE),
  #                        httr::progress())
  # disease_data <- httr::GET(com_dis,
  #                        httr::write_disk(disease_path, overwrite = TRUE),
  #                        httr::progress())
  # gene_content <- httr::content(gene_data, type = "text", encoding = "ISO-8859-1")
  # disease_content <- httr::content(disease_data, type = "text", encoding = "ISO-8859-1")
  compound_gene <- suppressWarnings(readr::read_tsv(com_gene,
                                                    show_col_types = FALSE,
                                                    name_repair = janitor::make_clean_names,
                                                    col_select = -c(number_input)
                                                    ))
  compound_disease <- suppressWarnings(readr::read_tsv(com_dis,
                                                       show_col_types = FALSE,
                                                       name_repair = janitor::make_clean_names,
                                                       col_select = -c(number_input)
                                                       ))
  cli_alert_success(c("Downloaded completed"))
  if (store) {
    output_dir <- destination(output_dir)
    gene_path <- glue::glue("{output_dir}/compound_gene.tsv")
    disease_path <- glue::glue("{output_dir}/compound_disease.tsv")
    compound_gene %>%
      readr::write_csv(file = gene_path)
    compound_disease %>%
      readr::write_csv(file = disease_path)
  }
  return(list(compound_gene = compound_gene, compound_disease = compound_disease))
}


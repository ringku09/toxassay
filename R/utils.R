
`%**%` <- function(x, y) matrixProd(x, y)

#' Function for matrix of one
#'
#' The function `mat_one()` is used to generate a matrix of one.
#'
#' @param n Number of rows
#' @param m Number of columns
#'
#' @return A matrix of ones5
#' @export
#'
#' @examples
#' mat_one(3, 4)
mat_one <- function(n, m) {
  matrix(1, nrow = n, ncol = m)
}


add_zero <- function(x) {
  if (nchar(x) == 1) {
    x <- paste0(0, x)
  } else x
  return(x)
}

#' Direct sum of two vectors or matrices
#'
#' The function ``
#' @param A A numeric matrix or vector
#' @param B A numeric matrix or vector
#'
#' @return
#' a matrix of direct sum
#' @export
#'
#' @examples
#' A <- matrix(1, nrow = 3, ncol = 3)
#' B <- matrix(2, nrow = 3, ncol = 3)
#' direct_sum(A, B)
direct_sum <- function(A, B) {
  if (is.matrix(A) & is.matrix(B)) {
    A <- A
    B <- B
  } else if (is.vector(A) & is.vector(B)) {
    A <- as.matrix(A)
    B <- as.matrix(B)
  } else stop( "Argument A and B must be a matrix or vector" )
  up_mat <- cbind(A, matrix(0, nrow = nrow(A), ncol = ncol(B)))
  down_mat <- cbind(matrix(0, nrow = nrow(B), ncol = ncol(A)), B)
  dir_sum <- rbind(up_mat, down_mat)
  return(dir_sum)
}

# # Set p-value sign based on coordinate
# plot_pval <- function(label, p_value) {
#   if (identical(label, "group")) {
#     p_val <- p_value
#   } else if (identical(label, "compound")) {
#     p_val <- p_value
#   } else if (identical(label, "dose")) {
#     p_val <- -p_value
#   } else if (identical(label, "time")) {
#     p_val <- -p_value
#   }
#   return(p_val)
# }

# # Set log FC sign based on coordinate
# plot_lfc <- function(label, logFC) {
#   if (identical(label, "group")) {
#     lfc <- logFC
#   } else if (identical(label, "compound")) {
#     lfc <- -logFC
#   } else if (identical(label, "dose")) {
#     lfc <- -logFC
#   } else if (identical(label, "time")) {
#     lfc <- logFC
#   }
#   return(lfc)
# }

# Find middle potion of bar plot for group to display text
gr_barpos <- function(gr_data, var) {
  zz <- gr_data %>%
    select({{var}})
  z <- as.vector(unlist(zz))
  if(all(z >= 0)) {
    txt_lev <- rev(cumsum(lag(rev(z), default = 0)) + rev(z)/2)
  } else if (any(z <  0) & any(z >= 0)) {
    txt_lev <- vector("numeric", length(z))
    z1 <- rev(cumsum(lag(rev(z[z >= 0]), default = 0)) + rev(z[z >= 0])/2)
    z2 <- cumsum(lag(rev(z[z < 0]), default = 0)) + rev(z[z < 0])/2
    txt_lev[z >= 0] <- z1
    txt_lev[z < 0] <- z2
  } else if (all(z < 0)) {
    txt_lev <- vector("numeric", length(z))
    z1 <- rev(cumsum(lag(rev(z[z >= 0]), default = 0)) + rev(z[z >= 0])/2)
    z2 <- cumsum(lag(rev(z[z < 0]), default = 0)) + rev(z[z < 0])/2
    txt_lev[z >= 0] <- z1
    txt_lev[z < 0] <- z2
    txt_lev <- rev(txt_lev)
  }
  return(txt_lev)
}

# Dataframe to matrix
df2matrix <- function(df, row_names = NULL) {
  df_mat <-  as.matrix(df)
  if (!is.null(rownames))
    rownames(df_mat) = row_names
  return(df_mat)
}

# Average fold change
avg_fc <- function(x, y, FC = TRUE) {
  avg_gr <- x %*% y$qAlfa
  if (FC) {
    fc <- avg_gr[1] - avg_gr[2]
  } else fc <- avg_gr[1]
  return(fc)
}

# destination path creator
destination <- function(output_dir) {
  if (is_missing(output_dir)) {
    output_dir <- tempdir()
  }
  if(!file.exists(output_dir)) {
    dir.create(output_dir)
  }
  return(output_dir)
}

#' @title Compute pairwise difference between matrix columns
#' @param x A data matrix of size n times p. Where rows are observations and
#' columns are features.
#' @return A matrix of size n times (p choose 2), where each column is the
#' difference between two of the original columns.
#' @importFrom magrittr set_colnames
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr imap
#' @export
#' @examples
#' n = 10
#' p = 3
#' x = matrix(rnorm(n*p), nrow = n, ncol = p, byrow = TRUE)
#' rownames(x) <- paste0("r", 1:n)
#' colnames(x) <- paste0("c", 1:p)
#' col_diff(x)
col_diff <- function(avg_data) {
  p <- ncol(avg_data)
  xin_list <- purrr::map(.x = seq_len(p - 1),
                         .f = ~ avg_data[, .x] - avg_data[, -c(seq_len(.x)), drop = FALSE])
  names(xin_list) <- colnames(avg_data)[-length(colnames(avg_data))]
  xin_list_name <- purrr::imap(
    .x = xin_list,
    .f = function(x, y) {
      mat_name = paste0(y, "-", colnames(x))
      colnames(x) = mat_name
      return(x)
    }
  )

  pair_diff <- do.call(cbind, xin_list_name)
  pair_mat <- tibble::as_tibble(pair_diff, rownames = "probe_id")
  return(pair_mat)
}


# with common legend
com_legend <- function(gg_plot) {
  tmp <- ggplot_gtable(ggplot_build(gg_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Capitalize first letter
block_fst <- function(x) {
  if (is.character(x)) {
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  } else x <- x
  return(x)
}


# compound name to abbreviation
comp_abbr <- function(comp_name, error_call = caller_env()) {
  if (is_empty(comp_name)) {
    cli_abort(c("{.var comp_name} must be non-empty.",
                "i" = "You have supplied an empty vector, please provide `name` or, `abbreviation` instred.")
              , call = error_call)
  }
  comp_tg <- tggates_com
  comp_is <- comp_name %in% comp_tg$compound_name
  if (all(comp_is)) {
    idx <- match(comp_name, comp_tg$compound_name)
    abbr <- comp_tg$compound_abbr[idx]
  }
  abbr_is <- comp_name %in% comp_tg$compound_abbr
  if (all(abbr_is)) {
    abbr <- comp_name
  }
  avail_com <- glue::englue("tggates_compounds()")
  if (any(comp_is) & any(abbr_is)) {
    cli_abort(c("Compound must be either in `name` or, `abbreviation`.",
                "x" = "Input {style_bold(col_blue(backtick(comp_name[comp_is])))} {?is/are} in compound name{?s} \\
                and {style_bold(col_green(backtick(comp_name [abbr_is])))} {?is/are} in compound abbreviation{?s}.",
                "i" = "Find available compound name/abbreviation in TG-GATEs by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
              , call = error_call)
  } else if (all(!comp_is) & all(!abbr_is)) {
    cli_abort(c("Compound must available in open TG-GATEs database.",
                "x" = "{style_bold(col_red(backtick(comp_name)))} compound{?s} name/abbrerviation \\
                {?is/are} not available in open TG-GATEs database.",
                "i" = "Please find the name of available compound by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
              , call = error_call)
  }else if (any(!comp_is) & all(!abbr_is)) {
    if (any(!comp_is)) {
      if (sum(comp_is) > length(comp_name)/2) {
        idx <- match(comp_name[comp_is], comp_tg$compound_name)
        abbr <- comp_tg$compound_abbr[idx]
        cli_alert_warning(c("Input {style_bold(col_red(backtick(comp_name[!comp_is])))} name{?s} of compound ",
                           "{?is/are} removed from the compound list, ",
                           "because {?this/these} input compound{?s} {?is/are} not available ",
                           "in open TG-GATEs database."),
                          wrap = TRUE)
      } else {
        cli_abort(c("Compound must available in open TG-GATEs database.",
                    "x" = "Input ({style_bold(col_red(backtick(comp_name[!comp_is])))}) name{?s} \\
                    of compound {?is/are} not available in open TG-GATEs database.",
                    "i" = "Please choose compound from the open TG-GATEs database,  \\
                    you can find compound list by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
                  , call = error_call)
      }
    }
  } else if (all(!comp_is) & any(!abbr_is)) {
    if (any(!abbr_is)) {
      if (sum(abbr_is) > length(comp_name)/2) {
        idx <- match(comp_name[abbr_is], comp_tg$compound_abbr)
        abbr <- comp_tg$compound_abbr[idx]
        cli_alert_warning(c("Input {style_bold(col_red(backtick(comp_name[!abbr_is])))} abbreviation{?s} of ",
                           "compound{?s} {?is/are} removed from the compound list, ",
                           "because {?this/these} compound{?s} {?is/are} not available ",
                           "in open TG-GATEs database."),
                          wrap = TRUE)
      } else {
        cli_abort(c("Compound must available in open TG-GATEs database.",
                    "x" = "Input ({style_bold(col_red(backtick(comp_name[!abbr_is])))}) \\
                    abbreviation{?s} of compound {?is/are} not available in open TG-GATEs database.",
                    "i" = "Please choose compound from the open TG-GATEs database,  \\
                    you can find available compound list by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
                  , call = error_call)
      }
    }
  }
  return(abbr)
}

# abbreviation or combination of abbreviation and compound name to compound name
abbr2name <- function(comp_abbr, error_call = caller_env()) {
  if (is_empty(comp_abbr)) {
    cli_abort(c("{.var comp_abbr} must be non-empty.",
                "i" = "You have supplied an empty vector, please provide `name` instred.")
              , call = error_call)
  }
  comp_tg <- compounds_tggates
  comp_out <- comp_abbr %in% c(comp_tg$compound_abbr, comp_tg$compound_name)
  avail_com <- englue("tggates_compounds()")
  if (any(!comp_out)) {
    cli_abort(c("Compound name must available in open TG-GATEs database.",
                "x" = "{style_bold(col_red(backtick(comp_abbr[!comp_out])))} compound{?s} name \\
                {?is/are} not available in open TG-GATEs database.",
                "i" = "Please find available compound name by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
              , call = error_call)
  }
  comp_is <- comp_abbr %in% comp_tg$compound_abbr
  if (all(comp_is)) {
    idx <- match(comp_abbr, comp_tg$compound_abbr)
    name <- comp_tg$compound_name[idx]
  }
  name_is <- comp_abbr %in% comp_tg$compound_name
  if (all(name_is)) {
    cli_alert_warning(c("Input {style_bold(col_red(backtick(comp_abbr)))} compound{?s}",
                        "{?is/are} already in name{?s},  no need to convert."),
                      wrap = TRUE)
    name <- comp_abbr
  } else if (any(name_is)) {
    cli_alert_warning(c("Input {.var comp_abbr} {style_bold(col_blue(backtick(comp_abbr[name_is])))} ",
                        "{?is/are} already in name{?s}"),
                      wrap = TRUE)
    rem_abbr <- comp_abbr[!name_is]
    rem_idx <- match(rem_abbr, comp_tg$compound_abbr)
    rem_name <- comp_tg$compound_name[rem_idx]
    name <- comp_abbr
    name[!name_is] <- rem_name
  }
  return(name)
}


# Find data type from attribute data
get_dtype <- function(attr_df, error_call = caller_env()) {
  species <- unique(attr_df$species)
  test_type <- unique(attr_df$test_type)
  organ <- unique(attr_df$organ_id)
  sin_rep <- unique(attr_df$sin_rep_type)
  if (length(species) > 1) {
    cli_abort(c("More than one `species` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} species {?is/are}  \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `species` to get appropiate data type.")
              , call = error_call)
  }
  if (length(test_type) > 1) {
    cli_abort(c("More than one `test_type` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} test type{?s} {?is/are}  \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `test_type` to get appropiate data type.")
              , call = error_call)
  }
  if (length(organ) > 1) {
    cli_abort(c("More than one `organ` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} organ{?s} {?is/are}  \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `organ` to get appropiate data type.")
              , call = error_call)
  }
  if (length(sin_rep) > 1) {
    cli_abort(c("More than one `sin_rep_type` not allowed.",
                "x" = "{style_bold(col_red(backtick(species)))} experiment type{?s} {?is/are} \\
                present in the attribute data.",
                "i" = "Please use attribute data with only one `sin_rep_type` to get appropiate data type.")
              , call = error_call)
  }
  dtype <- case_when(
    species == "Rat" & test_type == "in vivo" & organ == "Liver" & sin_rep == "Single" ~ "RVLS",
    species == "Rat" & test_type == "in vivo" & organ == "Liver" & sin_rep == "Repeat" ~ "RVLR",
    species == "Rat" & test_type == "in vivo" & organ == "Kidney" & sin_rep == "Single" ~ "RVKS",
    species == "Rat" & test_type == "in vivo" & organ == "Kidney" & sin_rep == "Repeat" ~ "RVKR",
    species == "Rat" & test_type == "in vitro" & organ == "Liver" & is.na(sin_rep) ~ "RVT",
    species == "Human" & test_type == "in vitro" & organ == "Liver" & is.na(sin_rep) ~ "HVT"
  )
  return(dtype)
}



# select_pcol <- function(pmat, lab = NULL, error_call = caller_env()) {
#     if (is.null(lab)) {
#       test_pmat <- pmat
#     } else {
#       lab <- gsub(" ", "", lab)
#       all_lab <- c("group",  "compound", "dose", "time")
#       unlist_lab <- unlist(strsplit(lab, "[^(A-Za-z)]"))
#       if (any(unlist_lab %in% "")) {
#         unlist_lab <- unlist_lab[-which(unlist_lab == "")]
#       }
#       if (!all(unlist_lab %in% all_lab)) {
#         miss_lab <- unlist_lab[!unlist_lab %in% all_lab]
#         remin_lab <- all_lab[!all_lab %in% unlist_lab]
#         cli_abort(c("The name of stage used in test level must be a valid stage.",
#                     "x" = "Input {style_bold(col_red(backtick(miss_lab)))} {?is/are} not a valid stage.",
#                     "i" = "Please use {style_bold(col_red(backtick(remin_lab)))} instread."),
#                   call = error_call)
#       }
#       join_sym <- unlist(strsplit(gsub("[^[:punct:]S]", "", lab), ""))
#       if (!all(join_sym %in% "&")) {
#         wrong_sym <- join_sym[!join_sym %in% "&"]
#         cli_abort(c("Test level condition contain wrong symbol.",
#                     "x" = "Input {style_bold(col_red(backtick(wrong_sym)))} {?is/are} not a valid symbol.",
#                     "i" = "Please use `&` instread."),
#                   call = error_call)
#       }
#       test_pmat <- as.matrix(pmat[, unlist_lab])
#     }
#     return(test_pmat)
#   }

# up-doen regulation calculator
up_down <- function(...,
                    probes,
                    ge_matrix,
                    metadata,
                    multicore = FALSE,
                    store = FALSE,
                    output_dir = missing_arg(),
                    error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  com_group <- test_group(...)
  output_dir <- destination(output_dir)
  compounds <- as.vector(unlist(com_group))
  ck_data <- update_data(
    compounds,
    ge_matrix = ge_matrix,
    metadata = metadata,
    probes = probes,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  ge_matrix <- ck_data$expression
  metadata <- ck_data$metadata
  lev_str <- data_str(comps_gr, metadata = metadata)
  design_mat <- get_matrix(lev_str)
  avg_gfc <- ge_matrix %**% design_mat$group_mat
  colnames(avg_gfc) <- names(comps_gr)
  rownames(avg_gfc) <- rownames(ge_matrix)
  up_downrg <- apply(avg_gfc, 1, function(x) ifelse(x[1]>x[2], "Upregulation", "Downregulation"))
  updown_df <- tibble::as_tibble(avg_gfc, rownames = "probe_id") %>%
    dplyr::mutate(Expression = factor(up_downrg,levels = c("Upregulation", "Downregulation")))
  return(updown_df)
}


split_two <- function(sentence) {
  words <- strsplit(sentence, "\\s+")[[1]]
  num_words <- length(words)

  if (num_words < 2) {
    first_half <- sentence
    second_half <- ""
  } else {
    mid_point <- ceiling(num_words / 2)
    first_half <- paste(words[1:mid_point], collapse = " ")
    second_half <- paste(words[(mid_point + 1):num_words], collapse = " ")
  }
  return(list(first_half = first_half, second_half = second_half))
}

split_sentence <- function(sentence, n_letters) {
  sentence <- as.character(sentence)
  words <- strsplit(sentence, " ")[[1]]
  word_lengths <- nchar(words)

  groups <- list()
  current_group <- c()
  current_sum <- 0

  for (i in seq_along(word_lengths)) {
    word_length <- word_lengths[i]
    if (current_sum + word_length > n_letters) {
      if (length(current_group) > 0) {
        groups[[length(groups) + 1]] <- paste(current_group, collapse = " ")
      }
      current_group <- c(words[i])
      current_sum <- word_length
    } else {
      current_group <- c(current_group, words[i])
      current_sum <- current_sum + word_length
    }
  }

  if (length(current_group) > 0) {
    groups[[length(groups) + 1]] <- paste(current_group, collapse = " ")
  }
  return(groups)
}



split_merge <- function(sentences, n_letters = 20) {
  split_x <- lapply(sentences, function(x) split_sentence(x, n_letters = n_letters))
  merge_split <- unlist(lapply(split_x, function(x) paste(unlist(x),collapse = "\n")))
  return(merge_split)
}

n_words <- function(sentence) {
  words <- strsplit(sentence, "\\s+")[[1]]
  n_word <- length(words)
  return(n_word)
}

# split_multiple <- function(sentence, nparts) {
#   words <- strsplit(sentence, "\\s+")[[1]]
#   n_char <- sapply(words, nchar)
#   cum_char <- cumsum(n_char)
#   for (i in 1:length(cum_char)) {
#     if(cum_char >= n_letters)
#   }
#
#
#   num_words <- length(words)
#
#   if (num_words < nparts) {
#     cli_abort(glue::glue("The sentence should contain at least {nparts} words."))
#   }
#
#   words_per_part <- ceiling(num_words / nparts)
#
#   sentence_parts <- vector("list", nparts)
#
#   for (i in 1:nparts) {
#     start_idx <- (i - 1) * words_per_part + 1
#     end_idx <- min(start_idx + words_per_part - 1, num_words)
#     sentence_parts[[i]] <- paste(words[start_idx:end_idx], collapse = " ")
#   }
#
#   return(sentence_parts)
# }


split2_merge <- function(sentences) {
  split_x <- lapply(sentences, function(x) split_two(x))
  max_char <- max(sapply(unlist(split_x), nchar))
  merge_split <- unlist(lapply(split_x,
                        function(x) {ifelse(x[[2]] != "",
                                            paste(unlist(x),collapse = "\n"), x[[1]])}))
  merge_split[nchar(sentences) < max_char] <- sentences[nchar(sentences) < max_char]
  return(merge_split)
}

title_tag <- function(title = NULL, title_size = rel(2), note = NULL,
                      note_size = rel(0.5), tag = NULL, tag_size = rel(1)) {
  if (!is.null(tag)) {
    text(x = -1,
         y = 1 ,
         labels = tag,
         cex = tag_size,
         font = 2)}
  if (!is.null(title)) {
    title(main = title,
          cex.main = title_size,
          font.main = 2
    )
  }
  if (!is.null(note)) {
    text(x = 1,
         y = -1 ,
         labels = note,
         pos = 2,
         cex = note_size
    ) #, family= "Times New Roman")
  }
}



rm_msg <- function(cat_msg) {
  cat("\r")
  cat(paste0(rep(" ", nchar(cat_msg)), collapse = ""))
  cat("\r")
}



get_class <- function(comps_gr, compounds) {
  class_assignment <- sapply(compounds, function(element) {
    matched_clusters <- sapply(comps_gr, function(cluster) element %in% cluster)
    match(TRUE, matched_clusters)
  })

  return(class_assignment)
}



round_up <- function(x, digits = 2) {
  multiplier <- 10^digits
  ifelse(x > 0, ceiling(x * multiplier) / multiplier, floor(x * multiplier) / multiplier)
}


add_or <- function(vec) {
  if (length(vec) == 1) {
    return(as.character(vec))
  }
  paste(paste(vec[-length(vec)], collapse = ", "), "or", vec[length(vec)])
}



check_internet <- function(site, time_out = 10) {
  if (!curl::has_internet()) {
    cli_abort(c("No internet connection.",
                "i" = "Please connect an internet and try again."))
  }
  response <- httr::GET(site, httr::timeout(time_out))
  if (!httr::status_code(response) == 200) {
    cli_abort(c("{gsub('_', ' ', deparse(substitute(site)))} FTP server not responding.",
                "i" = "Please try again later."))
  }
}


set_names <- function(level, metadata, name_column) {
  test_column(name_column, metadata)
  test_input(level, c("compound",
                      "compound-dose",
                      "compound-time",
                      "compound-dose-time"))
  lab <- metadata %>% dplyr::select({{name_column}}) %>% dplyr::pull()
  if (identical(level, "compound")) {
    col_names <- unique(lab)
  } else if (identical(level, "compound-dose")) {
    col_names <- unique(glue("{lab}.{metadata$dose_level}"))
  } else if (identical(level, "compound-time")) {
    col_names <- unique(glue("{lab}.{metadata$time_level}"))
  } else if (identical(level, "compound-dose-time")) {
    col_names <- unique(glue("{lab}.{metadata$dose_level}.{metadata$time_level}"))
  }
  return(col_names)
}

trans_color <- function(color, n_color=10) {
  lighten_color <- colorRampPalette(c(color, "#FFFFFF"))
  light_colors <- lighten_color(n_color+1)
  col_code <- light_colors[1:n_color]
  return(col_code)
}


color_light <- function(color, alpha) {
  lighten_color <- colorRampPalette(c(color, "#FFFFFF"))
  light_colors <- rev(lighten_color(11))
  col_code <- light_colors[alpha*10+1]
  return(col_code)
}

node_color <- function(n_nodes) {
  node_colors <- Polychrome::createPalette(n_nodes,
                                          c("#ff0000", "#00ff00", "#0000ff"))
  return(node_colors)
}

edge_color <- function(node_colors, alpha = 0.5) {
   edge_colors <- sapply(node_colors, FUN = color_light, alpha = alpha)
   return(edge_colors)
 }


#' Generate Random Samples from a Dirichlet Distribution
#'
#' This function generates random samples from a Dirichlet distribution using gamma function.
#'
#' @param n Integer. The number of samples to generate.
#' @param alpha Numeric vector. The concentration parameters for the Dirichlet distribution.
#'  Each element of `alpha` represents a parameter for the corresponding component of the Dirichlet distribution.
#'
#' @return A numeric vector of length `n * length(alpha)` representing random samples from the Dirichlet distribution.
#'
#' @examples
#' # Generate a random sample from a Dirichlet distribution with three components
#' rdirichlet(1, c(100, 100, 100))
#'
#' @export
rdirichlet <- function(n, alpha) {
  value <- rgamma(n * length(alpha), alpha, 1)
  rand_val <- value / sum(value)
  return(rand_val)
}

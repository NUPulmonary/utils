#' Clean GSEA Category Names
#'
#' Convert gene set enrichment analysis (GSEA) category names into readable
#' labels by replacing underscores with spaces, applying title case, preserving
#' specified acronyms, and correcting known spelling errors.
#'
#' @param pathway A character vector of category names, or a factor.
#'   Missing values are preserved.
#' @param drop_hallmark A single logical value indicating whether to remove
#'   the leading `"HALLMARK_"` prefix, ignoring case. Defaults to `TRUE`.
#' @param acronyms A character vector containing the desired capitalization of
#'   acronyms and other protected words. Matching is case-insensitive.
#'   Supply `character(0)` to disable acronym preservation.
#' @param corrections A named character vector of spelling corrections.
#'   Names are misspelled words; values are their replacements. Matching is
#'   case-insensitive and applies to whole words. Supply `character(0)` to
#'   disable spelling corrections.
#'
#' @details
#' Corrections are applied before capitalization. Words in `acronyms` retain
#' their specified capitalization; all other words are converted to title case.
#' Unknown acronyms cannot be inferred reliably from uppercase input and must
#' be added to the dictionary.
#'
#' The default dictionary covers acronyms and selected mixed-case labels in
#' the supplied Hallmark-style categories. Spelling correction is limited to
#' the explicit `corrections` dictionary; this function does not perform
#' general spell checking or validate gene set identifiers.
#'
#' Underscores become spaces, repeated whitespace is collapsed, and leading
#' and trailing whitespace is removed. Other punctuation is retained.
#' Direction labels are not expanded: `"UP"` becomes `"Up"` and `"DN"`
#' remains `"DN"` with the default dictionary.
#'
#' @return A character vector of the same length and order as `pathway`.
#'   Input names and missing values are preserved.
#'
#' @examples
#' pathways <- c(
#'   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#'   "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
#'   "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY",
#'   "Disease-Associated Microglia",
#'   "SenMayo",
#'   NA_character_
#' )
#'
#' clean_GSEA_names(pathways)
#' clean_GSEA_names(pathways, drop_hallmark = FALSE)
#'
#' # Supply a custom dictionary for another collection.
#' clean_GSEA_names(
#'   "CUSTOM_ATP_OXIGEN_RESPONSE",
#'   acronyms = c("ATP"),
#'   corrections = c(oxigen = "oxygen")
#' )
#'
#' @export
clean_GSEA_names <- function(
    pathway,
    drop_hallmark = TRUE,
    acronyms = c(
      "TNFA", "NFKB", "WNT", "TGF", "IL6", "JAK", "STAT3",
      "DNA", "G2M", "NOTCH", "PI3K", "AKT", "mTOR", "mTORC1",
      "E2F", "MYC", "V1", "V2", "p53", "UV", "IL2", "STAT5",
      "KRAS", "DN", "SenMayo"
    ),
    corrections = c(oxigen = "oxygen")
) {
  if (!is.character(pathway) && !is.factor(pathway)) {
    stop("`pathway` must be a character vector or factor.", call. = FALSE)
  }
  
  if (!is.logical(drop_hallmark) ||
      length(drop_hallmark) != 1L ||
      is.na(drop_hallmark)) {
    stop("`drop_hallmark` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.character(acronyms) || anyNA(acronyms)) {
    stop("`acronyms` must be a character vector without NA.", call. = FALSE)
  }
  
  if (!is.character(corrections) || anyNA(corrections)) {
    stop("`corrections` must be a character vector without NA.", call. = FALSE)
  }
  
  if (length(corrections) > 0L &&
      (is.null(names(corrections)) ||
       anyNA(names(corrections)) ||
       any(!nzchar(names(corrections))))) {
    stop("Every element of `corrections` must have a nonempty name.",
         call. = FALSE)
  }
  
  input_names <- names(pathway)
  x <- as.character(pathway)
  
  if (drop_hallmark) {
    x <- stringr::str_remove(
      x,
      stringr::regex("^HALLMARK_", ignore_case = TRUE)
    )
  }
  
  x <- stringr::str_squish(stringr::str_replace_all(x, "_", " "))
  
  correction_keys <- stringr::str_to_lower(names(corrections))
  acronym_keys <- stringr::str_to_lower(acronyms)
  
  # Process nonmissing elements explicitly to preserve NA values.
  present <- !is.na(x)
  x[present] <- stringr::str_replace_all(
    x[present],
    "\\b[[:alnum:]]+\\b",
    function(words) {
      correction_index <- match(
        stringr::str_to_lower(words), correction_keys
      )
      fix <- !is.na(correction_index)
      words[fix] <- unname(corrections[correction_index[fix]])
      
      result <- stringr::str_to_title(words)
      
      acronym_index <- match(
        stringr::str_to_lower(words), acronym_keys
      )
      keep <- !is.na(acronym_index)
      result[keep] <- unname(acronyms[acronym_index[keep]])
      
      result
    }
  )
  
  names(x) <- input_names
  x
}
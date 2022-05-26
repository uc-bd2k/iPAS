#' Produce a heatmap of iPAS pathway scores
#'
#' @param res output object from the iPAS_enrich function
#' @param row_fontsize font size for the row names (pathways), default is 6.
#' @param col_fontsize font size for the column names (input signatures),
#' default is 12.
#' @param short_names if TRUE, show shortened pathway names, not including
#' " - Homo sapiens (human)" at the end
#' @param pathway_names whether to show the names of KEGG pathways (e.g.
#' "mTOR signaling pathway - Homo sapiens (human)"), their unique identifiers
#' (e.g. "hsa04150") or both. "names" shows pathway names, "identifiers" shows
#' only ID numbers, "both" shows both, "none" shows no names.
#' @param category the categories of KEGG pathways for which to calculate iPAS
#' scores. We categorize pathways as "Disease", "Other", "Signaling", or
#' "Cancer". By default all categories/pathways are included.
#'
#'
#' @export

iPAS_heatmap = function(res, row_fontsize = 6, col_fontsize = 12, short_names = TRUE,
                       pathway_names = c("names", "identifiers", "both"),
                       category = c("Signaling", "Cancer", "Disease", "Other", "All")
                       ) {
  category = match.arg(category)
  pathway_names = match.arg(pathway_names)
  df = lapply(1:length(res), function(i) {
    name = names(res)[i]
    df = res[[i]]$ensemble
    df$sample = name
    df$null_dist = NULL
    return(df)
  }) %>% dplyr::bind_rows()

  if(short_names) {
    df$pathway_desc = df$pathway_desc %>% gsub(" - Homo sapiens \\(human\\)", "", .)
  }
  df$pathway_both = paste(df$pathway, df$pathway_desc)
  if(category != "All") {
    df = df %>% dplyr::filter(Category_New == category)
  }

  if(pathway_names == "names") {
    path_col = "pathway_desc"
  } else if(pathway_names == "identifiers") {
    path_col = "pathway"
  } else {
    path_col = "pathway_both"
  }
  f <- as.formula(paste(path_col, "~ sample"))

  mat = reshape2::acast(data = df, formula = f,
                        value.var = "z_score", fun.aggregate = mean)
  mat = mat[,names(res)]
  hm = ComplexHeatmap::Heatmap(mat, cluster_columns = T,
                               row_names_gp = grid::gpar(fontsize = row_fontsize),
                               column_names_gp = grid::gpar(fontsize = col_fontsize)
                               )
  return(hm)
}

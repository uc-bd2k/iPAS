#' Produce bar plots of top iPAS pathways
#'
#' This function creates bar plots that graph iPAS results. The input to the
#' function is the object returned from the iPAS_enrich function.
#'
#' @param res output object from the iPAS_enrich function
#' @param experiment the sample number or name of the sample from which to
#' graph results.
#' @param num_top the number of top pathways (by iPAS score) to plot
#' @param col_scale_max a positive number that defines the limits of the color
#' scale. The color scale is defined as c(-col_scale_max, col_scale_max).
#' @param pal scico package color palette. See here:
#' https://github.com/thomasp85/scico
#' @param pathway_names whether to show the names of KEGG pathways (e.g.
#' "mTOR signaling pathway - Homo sapiens (human)"), their unique identifiers
#' (e.g. "hsa04150") or both. "names" shows pathway names, "identifiers" shows
#' only ID numbers, "both" shows both, "none" shows no names.
#' @param short_names if TRUE, show shortened pathway names, not including
#' " - Homo sapiens (human)" at the end
#' @param rev T/F whether to reverse the order of the bars in the chart or not
#' @param fill the column to be used for the fill (color) of the bars. The
#' default is "z_score", but other options are "score" for the raw iPAS score
#' or "neg_log10_p" for the p-value (negative log10-transformed).
#' @param bar_height the column to be used for the height of the bars. The
#' choices are the same as above, the default is "neg_log10_p", the p-value
#' (negative log10-transformed).
#' @param arrange_col the column to use to order the bars, by default "p_emp"
#' for the empirical p-value (non-transformed).
#' @param category the categories of KEGG pathways for which to calculate iPAS
#' scores. We categorize pathways as "Disease", "Other", "Signaling", or
#' "Cancer". By default all categories/pathways are included.
#' @import ggplot2
#' @export

iPAS_bar = function(res, experiment = 1, num_top = 15, col_scale_max = 5,
                pal = "vikO",
                pathway_names = c("names", "identifiers", "both", "none"),
                short_names = TRUE,
				rev = FALSE,
                fill = c("z_score", "abs_z_score", "neg_log10_p"),
				bar_height = c("neg_log10_p", "z_score", "abs_z_score"),
				show_dir = c("Both", "Positive", "Negative"),
                arrange_col = c("z_score", "neg_log10_p", "p_emp", "abs_z_score"),
                category = c("Signaling", "Cancer", "Disease", "Other", "All")
) {
  category = match.arg(category, several.ok = TRUE)
  pathway_names = match.arg(pathway_names)
  fill = match.arg(fill)
  bar_height = match.arg(bar_height)
  arrange_col = match.arg(arrange_col)
  show_dir = match.arg(show_dir)
  if(is.numeric(experiment)) {
    exper = names(res)[experiment]
  } else {
    exper = experiment
  }
  ens_df = res[[experiment]]$ensemble
  if(!identical(category, "All")) {
    ens_df = ens_df %>% dplyr::filter(Category_New %in% category)
  }
  if(short_names) {
    ens_df$pathway_desc = ens_df$pathway_desc %>% gsub(" - Homo sapiens \\(human\\)", "", .)
  }
  ens_df$pathway_both = paste(ens_df$pathway, ens_df$pathway_desc)
  df = ens_df %>%
    dplyr::mutate(neg_log10_p = -log10(p_emp)) %>%
  	dplyr::mutate(abs_z_score = abs(z_score))# %>%
    ### select most significant pathways
    #dplyr::arrange( p_emp ) %>%
    #dplyr::slice_head(n = num_top)
  
  if(show_dir == "Both") {
  	show_dir = c("Positive", "Negative")
  }
  df = df %>% dplyr::filter(direction %in% show_dir)
  
  if(arrange_col %in% c("neg_log10_p", "z_score", "abs_z_score")) {
  	df = df %>% dplyr::arrange( desc(!!sym(arrange_col)) )
  } else {
  	### for arranging by empirical p-value
  	df = df %>% dplyr::arrange( !!sym(arrange_col) )
  }
  if(rev) { df = df[dim(df)[1]:1,] }
  
  ### look at only top n pathways
  df = df %>% dplyr::slice_head(n = num_top)

  if(pathway_names == "names") {
    path_col = "pathway_desc"
  } else if(pathway_names == "identifiers") {
    path_col = "pathway"
  } else if(pathway_names == "both") {
    path_col = "pathway_both"
  } else {
    path_col = "pathway"
  }

  df[[path_col]] = factor(df[[path_col]], levels = rev(as.character(df[[path_col]])))
  #print(head(levels(df[[path_col]])))
  
  xmax = max(df[[bar_height]], na.rm = T)
  xmin = min(df[[bar_height]], na.rm = T)
  ### make xmax and xmin at most/least 3 and -3
  xmax = max(xmax, 3)
  xmin = min(xmin, -3)
  col_max = max(abs(df[[fill]]), na.rm = T)
  ### make col_max at least 3
  col_max = max(col_max, 3)
  if(bar_height == "neg_log10_p") {
  	xlims = c(0, xmax)
  } else {
  	xlims = c(xmin, xmax)
  }
  if(fill == "neg_log10_p") {
  	col_lims = c(0, col_max)
  } else {
  	col_lims = c(-col_max, col_max)
  }
  # gg_base = ggplot(df) + 
  # 	geom_col(aes(x = !!sym(path_col), 
  # 				 y = !!sym(bar_height), 
  # 				 fill = !!sym(fill))
  # 			 )
  gg_base = ggplot(df) + 
  	geom_segment(aes(x = !!sym(path_col),
  					 xend = !!sym(path_col),
  				 y = !!sym(bar_height),
  				 yend = 0,
  				 colour = !!sym(fill)
  				 ),
  				 size = 7
  	)
  if(!pathway_names == "none") {
    gg1 = gg_base
  } else {
    gg1 = gg_base +
      theme(axis.text.y = element_blank(), legend.position = "none") 
  }
  gg2 = gg1 +
  	#scico::scale_fill_scico(palette = pal, limit = col_lims)
  	scico::scale_color_scico(palette = pal, limit = col_lims) +
  	coord_flip(ylim = xlims) +
  	theme_bw() +
  	labs(title = exper)
  if(bar_height == "neg_log10_p") {
    gg2 = gg2 + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  }
  return(gg2)
}

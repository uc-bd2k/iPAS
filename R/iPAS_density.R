#' Produce a density plot of iPAS permutation scores for a given pathway
#'
#' This function returns a plot of the probability density of the null
#' distribution of iPAS scores for permutation signatures for a given pathway.
#' The iPAS score for the (non-permuted) input signature is shown by a red bar.
#'
#' @param res output object from the iPAS_enrich function
#' @param experiment the sample number or name of the sample from which to
#' graph results.
#' @param score_type the type of iPAS score to graph. "score" will show raw iPAS
#' scores and "z_score" will show normalized iPAS scores.
#' @param path the KEGG pathway ID of the pathway to graph (e.g. "hsa04150" for
#' the mTOR pathway)
#' @param fill the color to fill in the density plot, "light blue" by default.
#'
#' @export

iPAS_density = function(res, experiment, score_type = c("z_score", "score"),
                        path, fill = "light blue") {
  score_type = match.arg(score_type)
  if(is.numeric(experiment)) {
    exper = names(res)[experiment]
  } else {
    exper = experiment
  }
  ens_df = res[[experiment]]$ensemble
  path_df = ens_df %>% dplyr::filter(pathway == path)

  if(score_type == "z_score") {
    null_dist = path_df$null_dist[[1]]
    null_dist = (null_dist - mean(null_dist, na.rm = T))/sd(null_dist, na.rm = T)
    z = path_df$z_score[1]
  } else if(score_type == "score") {
    null_dist = path_df$null_dist[[1]]
    z = path_df$score[1]
  } else {
    stop("score_type is not available")
  }
  df = data.frame(pathway = path, null_dist = null_dist)

  p = path_df$p_emp[1]
  ann = paste0("p-value: ", p)
  annotations <- data.frame(
    xpos = c(-Inf,-Inf,Inf,Inf),
    ypos =  c(-Inf, Inf,-Inf,Inf),
    annotateText = c("","","",ann),
    hjustvar = c(0,0,1,2) ,
    vjustvar = c(0,1.0,0,2))

  gg = ggplot(df) + geom_density(aes(x=null_dist), fill = fill) +
    labs(title = exper, subtitle = path) + geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                                            vjust=vjustvar,label=annotateText)) +
    geom_vline(xintercept = z, colour = "red") + theme_classic()

  return(gg)
}

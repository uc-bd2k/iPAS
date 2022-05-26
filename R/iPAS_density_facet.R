#' Produce density plots of iPAS permutation scores for a given pathway for
#' multiple experiments
#'
#' This function is the same as the iPAS_density function, except for all
#' experiments at once.
#' The function returns a plot of the probability density of the null
#' distribution of iPAS scores for permutation signatures for a given pathway.
#' The iPAS score for the (non-permuted) input signature is shown by a red bar.
#'
#' @param res output object from the iPAS_enrich function
#' @param score_type the type of iPAS score to graph. "score" will show raw iPAS
#' scores and "z_score" will show normalized iPAS scores.
#' @param path the KEGG pathway ID of the pathway to graph (e.g. "hsa04150" for
#' the mTOR pathway)
#' @param fill the color to fill in the density plot, "light blue" by default.
#'
#' @import ggplot2
#' @export

iPAS_density_facet = function(res, score_type = c("z_score", "score"), path,
                              fill = "light blue") {
  score_type = match.arg(score_type)
  expers = names(res)

  null_df = lapply(1:length(res), function(i) {
    name = names(res)[i]
    df = res[[i]]$ensemble
    df$sample = name

    path_df = df %>% dplyr::filter(pathway == path)

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

    sample = path_df$sample[1]
    p = path_df$p_emp[1]
    df1 = data.frame(pathway = path, null_dist = null_dist, sample = sample, p = p, z = z)
    return(df1)
  }) %>% dplyr::bind_rows()

  null_df$sample = factor(null_df$sample, levels = names(res))

  vline_df =  null_df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(z = unique(z), p = unique(p))
  vline_df = vline_df[match(names(res), as.character(vline_df$sample)),]
  vline_df$sample = factor(vline_df$sample, levels = names(res))

  vline_df$ann = paste0("p-value: ", vline_df$p %>% prettyNum(digits = 2))

  gg = ggplot(null_df) + geom_density(aes(x=null_dist), fill = fill) +
    labs(title = path) +
    # geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
    #           vjust=vjustvar,label=annotateText)) +
    facet_grid(sample ~ .) + theme_classic() +
    geom_vline(data=vline_df, aes(xintercept=z), colour="red")

  gg = egg::tag_facet(gg,
            x = Inf, y = Inf,
            vjust = 1, hjust = 1.5,
            open = "", close = "",
            fontface = 1,
            size = 4,
            #family = "arial",
            tag_pool = vline_df$ann)

  gg = egg::tag_facet(gg,
                      x = -Inf, y = Inf,
                      vjust = 1, hjust = -0.1,
                      open = "", close = "",
                      fontface = 1,
                      size = 4,
                      #family = "arial",
                      tag_pool = vline_df$sample)

  return(gg)
}

#' Calculate iPAS pathway scores for a given input signature
#'
#' This function takes in a query signature (ex: log fold-change values of gene
#' expression for a perturbation experiment) and calculates iPAS scores for
#' each KEGG pathway. Many query signatures may be submitted at once. The input
#' signature(s) must be in the form of a matrix where each column is a signature
#' and each row is a gene. Row names must be Entrez IDs. Names for each
#' signature may be given as column names.
#'
#' @param query a matrix with the input/query signatures to perform pathway analysis on.
#' The row names should be genes (Entrez IDs) and the column names should be
#' names for the samples/signatures.
#' @param similarity the type of similarity measure to use for iPAS. Choices are
#' "Pearson" for Pearson correlation, "cosine" for cosine similarity, and
#' "dot_product" for the dot product.
#' @param gene_type the type of names used for the genes. Currently only Entrez
#' IDs are available.
#' @param category the categories of KEGG pathways for which to calculate iPAS
#' scores. We categorize pathways as "Disease", "Other", "Signaling", or
#' "Cancer". By default all categories/pathways are included.
#' @param perm the number of permutations to perform (1000 by default)
#' @param testing T/F value or integer. If TRUE, only calculate iPAS scores for
#' the first 5 pathways. If an integer n, calculate iPAS scores for the first n
#' pathways. Used in testing.
#' @param return_individual_cl T/F value, whether to return similarity scores
#' for each individual cell line
#' @param return_null_dist T/F value, whether to return the null distribution
#' of iPAS scores for each pathway (permutation scores). TRUE is required for
#' the iPAS_density and iPAS_density_facet functions, which graph the density
#' of the null distribution. FALSE by default.
#' @param overlap_min the minimum number of overlapping genes required between
#' the input signature and a pathway signature to calculate similarity (default
#' is 10).
#' @param ncores number of cores to use for calculation (via parallel package),
#' default is 1.
#' @param seed a seed (integer value) to use for random number generation,
#' passed to set.seed. This can be used to make the analysis reproducible, since
#' it involved random permutation.
#' @param print_updates T/F value, whether to print updates while calculating,
#' FALSE by default.
#'
#' @return a data frame (tibble) with correlation and p-value results for the input signature
#' compared with each PAS signature
#' @export

iPAS_enrich = function(query,
                        similarity = c("Pearson", "cosine", "dot_product"),
                        gene_type = "entrez",
                        category = c("Disease", "Other", "Signaling", "Cancer"),
                        perm = 1000,
                        testing = F,
                        return_individual_cl = T,
                        return_null_dist = F,
                        overlap_min = 10,
                        ncores = 1,
                        seed = NULL,
                        print_updates = F,
                        updateProgress = NULL
                        ) {

  call = list(query = query, similarity = similarity, gene_type = gene_type,
              category = category, perm = perm, testing = testing,
              return_individual_cl = return_individual_cl,
              return_null_dist = return_null_dist,
              ncores = ncores,
              seed = seed
              )
# Check inputs ------------------------------------------------------------
  similarity = match.arg(similarity)
  gene_type = match.arg(gene_type)
  #cell_line = match.arg(cell_line, several.ok = TRUE)
  category = match.arg(category, several.ok = TRUE)
  all_cl = c("A375", "A549", "HA1E", "HCC515",
             "HEKTE", "HEPG2", "HT29", "MCF7",
             "NPC", "PC3", "SW480", "VCAP")

# Load pathway metadata and PAS -------------------------------------------

  #### Load pathway metadata
  path_meta = iPAS::kegg_meta %>% dplyr::rename(pathway = Pathway)
  ens_weights = iPAS::ens_weights_no_intercept
  ens_pathways = colnames(ens_weights)

  ## Filter to only the specified categories
  path_meta = path_meta %>% dplyr::filter(Category_New %in% category)

  all_paths = path_meta$pathway
  all_paths = intersect(all_paths, ens_pathways)
  if(testing) {
    if(is.logical(testing)) {
      all_paths = all_paths[1:5]
    }
    if(is.numeric(testing)) {
      all_paths = all_paths[1:testing]
    }
  }

  #### Load ensemble signatures
  pas_by_path = iPAS::iPAS_by_path


  l1k_genes = pas_by_path[[1]] %>% rownames()

  #pas_ls = iPAS::pas_ls
  #ens_weights = readRDS("/opt/raid10/genomics/nick/Projects/Git/MTL/output/ensemble_beta_w2_0.rds")
  #ens_weights = ens_weights[-1,]
  ##check1 = all.equal(rownames(ens_weights), names(pas_ls))

# Loop over each sample of input/query ------------------------------------
  res = parallel::mclapply(1:dim(query)[2], function(i) {
    if(print_updates) {
      print(i)
    }
    tmp_query = query[,i]
    tmp_query = tmp_query[!is.na(tmp_query)]
    tmp_query = tmp_query[!is.infinite(tmp_query)]
    ### Create null distribution signatures
    if(!is.null(seed)) set.seed(seed)
    if(length(tmp_query) >= 978) {
      query_perm_full = sapply(1:perm, function(i) {
        sample(tmp_query, size = 978, replace = FALSE)
      })
    } else {
      query_perm_full = sapply(1:perm, function(i) {
        sample(tmp_query, size = 978, replace = TRUE)
      })
    }

  # Loop over each each pathway ---------------------------------------------
    #tictoc::tic()
    tmp_res = lapply(1:length(all_paths), function(j){
      x = all_paths[j]
      if(print_updates) {
        print(x)
      }
      # If we were passed a progress update function, call it
      if (is.function(updateProgress)) {
        n_samp = dim(query)[2]
        n_path = length(all_paths)
        total = n_samp*n_path
        prog = (i - 1)*n_path + j
        val = prog/total
        text <- paste0("Sample ", i, " of ", n_samp," - ", "Pathway ", j, " of ", n_path)
        updateProgress(value = val, detail = text)
      }
      pas_tmp = pas_by_path[[x]]
      null_check = is.null(pas_tmp)
      if(null_check) {
        print("NULL check")
        return(NULL)
      }
      ### set rownames of permutations to L1000 genes
      rownames(query_perm_full) = rownames(pas_tmp)
      check_tmp = all.equal(colnames(pas_tmp), all_cl)
      if(!check_tmp) stop("columns of 'pas_tmp' matrix out of order: iPAS_enrich function")
      # Loop over each cell line ------------------------------------------------
      tmp_cl = lapply(all_cl, function(cl) {
        if(print_updates) {
          #print(cl)
        }
        pas_ind = pas_tmp[,cl]
        ### match genes from tmp_query with those from pas_tmp and remove NAs
        pas_ind = pas_ind[!is.na(pas_ind)]
        pas_ind = pas_ind[!is.infinite(pas_ind)]

        shared_genes = intersect(names(pas_ind), names(tmp_query))
        overlap = length(shared_genes)
        pas_ind = pas_ind[shared_genes]
        query_ind = tmp_query[shared_genes]
        query_perm = query_perm_full[shared_genes, , drop = F]
        ### send to pathway scoring function
        wt = ens_weights[cl,x]
        sim_tmp = iPAS_score_perm(query = query_ind, query_perm = query_perm, pas = pas_ind, similarity = similarity, perm = perm, pathway = x, cell_line = cl, ens_weight = wt, overlap_min = overlap_min)
        return(sim_tmp)
      }) %>% magrittr::set_names(all_cl)

      # Calculate ensemble pathway scores ---------------------------------------
      score_vec = sapply(tmp_cl, function(x) x$score)
      ens_tmp = ens_weights[,x]
      score = sum(score_vec*ens_tmp)

      gene_scores_ind = sapply(tmp_cl, function(x) {
        tmp = x$gene_scores
        gene_vec = rep(0, 978) %>% magrittr::set_names(l1k_genes)
        gene_vec[match(names(tmp), names(gene_vec))] = tmp
        return(gene_vec)
      })

      gene_scores_ens_mat = gene_scores_ind %*% as.matrix(ens_tmp)
      gene_scores_ens = gene_scores_ens_mat[,1]

      check_tmp = all.equal(names(score_vec), names(ens_tmp))
      if(!check_tmp) stop("cell line names not matching")

      perm_scores_cl = sapply(tmp_cl, function(x) x$null_dist)
      perm_scores = perm_scores_cl %*% as.matrix(ens_tmp, ncol = 1)

      ecdf_perm = ecdf(perm_scores)
      mean_perm = mean(perm_scores, na.rm = TRUE)
      sd_perm = sd(perm_scores, na.rm = TRUE)
      z = (score - mean_perm)/sd_perm
      if(sign(score) > 0) {
        direction = "Positive"
        p_emp = (sum(perm_scores >= score) + 1)/perm ## add 1 to prevent p = 0
      } else if (sign(score) < 0) {
        direction = "Negative"
        p_emp = (sum(perm_scores <= score) + 1)/perm ## add 1 to prevent p = 0
      } else {
        direction = "Zero"
        p_emp = 1
      }
      overlap = lapply(tmp_cl, function(x) x$overlap)

      out_list = list(
        ensemble = list(score = score, z_score = z, p_emp = p_emp,
          overlap = overlap, direction = direction, similarity = similarity,
          gene_scores_ens = gene_scores_ens,
          pas_ens = list(pas_tmp),
          query_l1k = list(tmp_query),
          pathway = x, cell_line = NA, null_ecdf = ecdf_perm,
          null_dist = perm_scores, null_mean = mean_perm, null_sd = sd_perm
          )
      )
      ens_and_cl = c(tmp_cl, out_list)
      #### format list as a tibble and return
      out_tbl_ens = tibble::tibble(pathway = x, cell_line = "ensemble",
                                   score = round(score, 3),
                               z_score = round(z, 3), p_emp = p_emp, direction = direction,
                               ens_weight = NA, gene_scores_ens = list(gene_scores_ens),
                               score_type = similarity, overlap = list(overlap), null_dist = list(perm_scores)
      )
      ind_tbls = lapply(tmp_cl, function(x) x$tbl) %>% dplyr::bind_rows()
      full_tbl = ind_tbls %>% dplyr::bind_rows(out_tbl_ens)

      return( list(tbl = full_tbl) )
    }) %>% magrittr::set_names(all_paths)
    #tictoc::toc()
    tbl = lapply(tmp_res, function(x) x$tbl) %>% dplyr::bind_rows()
    if(!return_null_dist) tbl$null_dist = NULL
    ### merge table with pathway metadata
    tbl2 = tbl %>%
      dplyr::filter(cell_line == "ensemble") %>%
      dplyr::left_join(path_meta, by = "pathway") %>%
      dplyr::arrange( p_emp, desc(abs(z_score))) %>%
      dplyr::relocate(pathway, pathway_desc, score, z_score, p_emp,
                      direction, cell_line, Category_Lv1,
                      Category_Lv2, Category_New, size, path_nodes_list,
                      overlap
                      )

    if(return_individual_cl) {
      pathway_ind = tbl %>%
        dplyr::filter(cell_line != "ensemble") %>%
        split(.$pathway)
    } else {
      pathway_ind = tibble::tibble()
    }
    res_ls = list(ensemble = tbl2, cl = pathway_ind)
    return(res_ls)
  }, mc.cores = ncores) %>% magrittr::set_names(colnames(query))

  return(res)
}

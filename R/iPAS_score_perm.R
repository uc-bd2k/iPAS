#' Create a null distribution of iPAS pathway scores using permutations
#'
#'

iPAS_score_perm = function(query, query_perm, pas, similarity, perm, pathway, cell_line, ens_weight, overlap_min) {

  call = list(query = query, query_perm = query_perm, pas = pas, similarity = similarity, perm = perm, pathway = pathway, cell_line = cell_line, ens_weight = ens_weight)
  N = length(query)
  if(N < overlap_min) {
    out_tbl = tibble::tibble(pathway = pathway, cell_line = cell_line, score = 0,
                             z_score = 0, p_emp = 1, direction = "Zero", ens_weight = round(ens_weight, 3),
                             score_type = similarity, overlap = list(0), null_dist = list(rep(0, perm))
    )

    out_list = list(score = 0, z_score = 0, p_emp = 1, overlap = 0,
                    direction = "Neither", similarity = similarity,
                    pathway = pathway, cell_line = cell_line,
                    null_ecdf = ecdf(rep(0, perm)), null_dist = rep(0, perm),
                    null_mean = 0, null_sd = NA, call = call, tbl = out_tbl)
    return(out_list)
  }

  ### Score query signature and null distribution signatures
  if(similarity == "Pearson") {
    z1 = scale(query, center = TRUE)
    z2 = scale(pas, center = TRUE)

    z3 = scale(query_perm, center = TRUE)
    #score_test = (1/(N-1))*sum(z1*z2)
    #perm_scores_test = (1/(N-1))*as.vector( t(z3) %*% as.matrix(z2) )
    #perm_scores_test = (1/(N-1))*( colSums( apply(z3, 2, function(x) {x*z2}) ) )

    score = cor(query, pas)
    perm_scores = cor(query_perm, pas) %>% as.vector()

    gene_scores = (1/(N-1))*(z1*z2)
    gene_scores_perm = (1/(N-1))*( apply(z3, 2, function(x) {x*z2}) )
  } else if(similarity == "cosine") {
    z1 = scale(query, center = FALSE)
    z2 = scale(pas, center = FALSE)

    z3 = scale(query_perm, center = TRUE)

    #z1 = scale(query_perm[,1], center = FALSE) ### line for testing
    score = (1/(N-1))*sum(z1*z2)
    perm_scores = t(scale(query_perm, center = FALSE)) %*% as.matrix(z2, ncol = 1)
    perm_scores = (1/(N-1))*as.vector(perm_scores)

    gene_scores = (1/(N-1))*(z1*z2)
    gene_scores_perm = (1/(N-1))*( apply(z3, 2, function(x) {x*z2}) )
  } else if(similarity == "dot_product") {
    z1 = query
    z2 = pas
    #z1 = query_perm[,1] ### line for testing
    score = (1/(N-1))*sum(z1*z2)
    perm_scores = t(query_perm) %*% as.matrix(z2, ncol = 1)
    perm_scores = (1/(N-1))*as.vector(perm_scores)

    gene_scores = (1/(N-1))*(z1*z2)
    gene_scores_perm = (1/(N-1))*( apply(z3, 2, function(x) {x*z2}) )
  }

  gene_scores = gene_scores[,1]
  ### testing alternative Pearson calculation
  # if(similarity == "Pearson2") {
  #   z1 = scale(query)
  #   z2 = scale(pas)
  #   N = length(z1)
  #   score = (1/(N-1))*sum(z1*z2)
  # }
  #print(perm_scores)
  ecdf_perm = ecdf(perm_scores)
  mean_perm = mean(perm_scores, na.rm = TRUE)
  sd_perm = sd(perm_scores, na.rm = TRUE)
  z = (score - mean_perm)/sd_perm
  if(sign(score) > 0) {
    direction = "Positive"
    p_emp = sum(perm_scores >= score)/perm
  } else if (sign(score) < 0) {
    direction = "Negative"
    p_emp = sum(perm_scores <= score)/perm
  } else {
    direction = "Neither"
    p_emp = 1
  }

  out_tbl = tibble::tibble(pathway = pathway, cell_line = cell_line, score = round(score, 3),
                           z_score = round(z, 3), p_emp = p_emp, direction = direction,
                           gene_scores = list(gene_scores),
                           query_scaled = list(z1),
                           pas_scaled = list(z2),
                           ens_weight = round(ens_weight, 3),
                           score_type = similarity, overlap = list(N), null_dist = list(perm_scores)
                           )

  out_list = list(score = score, z_score = z, p_emp = p_emp, overlap = N,
                  direction = direction, similarity = similarity,
                  gene_scores = gene_scores,
                  query_scaled = z1,
                  pas_scaled = z2,
                  pathway = pathway, cell_line = cell_line,
                  null_ecdf = ecdf_perm, null_dist = perm_scores,
                  null_mean = mean_perm, null_sd = sd_perm, call = call, tbl = out_tbl)
  return(out_list)
}

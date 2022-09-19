
library(tidyverse)


gene_set_enrichment = function(candidate_gene_sets, reference_gene_sets, all_genes=NULL) {
    if (is.null(all_genes)) {
        all_genes = sort(unique(unlist(reference_gene_sets)))
    }
    candidates = gene_set_to_matrix(candidate_gene_sets, all_genes)
    positives = gene_set_to_matrix(reference_gene_sets, all_genes)
    result = hyper_test(candidates, positives)
    return(result)
}

gene_set_to_matrix = function(gene_set, all_genes = sort(unique(unlist(gene_set)))) {
    go_sets = enframe(gene_set, name = "term", value = "gene") %>%
        unnest("gene") %>%
        filter(gene %in% all_genes)

    gene_symbols = factor(go_sets$gene, levels = all_genes)
    go_terms = as.factor(go_sets$term)

    result = Matrix::sparseMatrix(as.numeric(gene_symbols), as.numeric(go_terms),
                            dim = c(nlevels(gene_symbols), nlevels(go_terms)),
                          dimnames = list(levels(gene_symbols), levels(go_terms)))
    return(1*result)
}

hyper_test = function(candidates, positives) {
    white_drawn = as.matrix(Matrix::crossprod(positives, candidates))
    white_total = Matrix::colSums(positives)
    black_total = nrow(positives) - white_total
    drawn_total = matrix(Matrix::colSums(candidates), ncol(positives), ncol(candidates), byrow = TRUE)
    # phyper() computes the proportion of scenarios that are *strictly* worse,
    # the -0.5 correction makes sure that the equally worse scenario is included
    p_value = phyper(white_drawn - 0.5, white_total, black_total, drawn_total, lower.tail = FALSE)
    fdr = matrix(p.adjust(p_value, "fdr"), nrow = nrow(p_value), dimnames = dimnames(p_value))
    #fold_enrichment = white_drawn * (white_total + black_total) / (white_total * drawn_total)
    white_not_drawn = white_total - white_drawn
    black_drawn = drawn_total - white_drawn
    black_not_drawn = black_total - black_drawn
    odds_ratio = white_drawn * black_not_drawn / (black_drawn * white_not_drawn)
    result = data.frame(
        candidate_set = rep(colnames(white_drawn), each=nrow(white_drawn)),
        candidate_size = rep(Matrix::colSums(candidates), each=nrow(white_drawn)),
        reference_set = rep(rownames(white_drawn), ncol(white_drawn)),
        reference_size = rep(white_total, ncol(white_drawn)),
        overlap = as.vector(white_drawn),
        odds_ratio = as.vector(odds_ratio),
        p_value = as.vector(p_value),
        fdr = as.vector(fdr)
    )
    return(result)
}


#' predictors is a matrix where each column is a predictor and each row is a sample.
#' label_matrix is a binary matrix where columns are labels and each row is a sample.
#' 1 indicates that the sample on this row belongs to the label on this column.
compute_aurocs = function(predictors, label_matrix, return_tie_correction = FALSE) {
    label_matrix = as.matrix(label_matrix)
    n_positives = colSums(label_matrix)
    n_negatives = nrow(label_matrix) - n_positives
    if (is(predictors, "dgCMatrix")) {
        # we shift all ranks after the matrix multiplication to keep
        # the predictor matrix sparse
        ranks = rank_sparse(predictors)
        sum_of_positive_ranks = as.matrix(crossprod(label_matrix, ranks)) +
            outer(n_positives, rank_zero(predictors))
    } else {
        predictors = as.matrix(predictors)
        ranks = matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
        sum_of_positive_ranks = crossprod(label_matrix, ranks)
        colnames(sum_of_positive_ranks) = colnames(predictors)
    }
    if (return_tie_correction) {
      tie_correction = compute_tie_correction(ranks)
    }
    result = (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
    if (return_tie_correction) {
      return(list(aurocs = result, tie_corrections = tie_correction))
    } else {
      return(result)
    }
}

design_matrix = function(cell_type, scale_columns=FALSE) {
  factors = levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result = model.matrix(~as.factor(cell_type)-1)
  } else {
    result = matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) = factors
  if (scale_columns) {
    result = scale(result, center = FALSE, scale = colSums(result))
  }
  return(result)
}

# these 2 functions only rank non-zeros, implicitly shifting the matrix of ranks
# to keep the matrix sparse according to the formula:
#   true_ranks(M) = rank_zero(M) + rank_sparse(M), where:
#     + rank_zero(M) = (#zeros + 1)/2 = ((nrow(M) - diff(M@p)) + 1) / 2
#     + rank_sparse(M) = rank(nonzero element) + (#zeros - 1)/2
# faster than solution using base::tapply, probably data.table would be faster
rank_sparse = function(M) {
    nnz = diff(M@p)
    ranks = tibble(x = M@x, j = rep.int(1:ncol(M), nnz)) %>%
        group_by(j) %>%
        mutate(rank_ = colRanks(as.matrix(x), ties.method = "average"))
    R = M
    R@x = ranks$rank_ + rep.int((nrow(M)-nnz-1)/2, nnz)
    return(R)
}

rank_zero = function(M) {
    return(((nrow(M) - diff(M@p)) + 1) / 2)
}

# For the following two functions, see
#
#   https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
#
# The tie correction effectively computes lost variance because of ties (compared to discrete uniform).
# Computing the Wikipedia formula naively is slow, this method is equivalent and fast.
compute_tie_correction = function(ranks) {
    if (is(ranks, "dgCMatrix")) {
        observed_var = colVars_sparse(ranks)
    } else {
        observed_var = matrixStats::colVars(as.matrix(ranks))
    }
    max_var = var(seq_len(nrow(ranks)))
    return((max_var-observed_var) * 12 / nrow(ranks))
}

auroc_p_value = function(aurocs, label_matrix, two_tailed = TRUE, tie_correction = 0, log.p = FALSE) {
    p = colSums(label_matrix)
    n = nrow(label_matrix) - p
  
    # Careful: NAs arise from tie_correction (predictor with 0 variance)
    if (length(tie_correction) > 1) {
        Z = (aurocs - 0.5) * sqrt(12*n*p)
        Z = t(t(Z) / sqrt(nrow(label_matrix)+1-tie_correction))
    } else {
        Z = (aurocs - 0.5) / sqrt((nrow(label_matrix)+1-tie_correction)/(12*n*p))
    }
  
    result = Z
    if (two_tailed) {
        is_not_na = !is.na(Z)
        result[Z<=0 & is_not_na] = pnorm(Z[Z<=0 & is_not_na], log.p = log.p) * 2
        result[Z>0 & is_not_na] = pnorm(Z[Z>0 & is_not_na], lower.tail=FALSE, log.p = log.p) * 2
    } else {
        result = pnorm(Z, lower.tail=FALSE, log.p = log.p)
    }
    return(result)
}

colVars_sparse = function(M) {
    result = (Matrix::colMeans(M**2) - Matrix::colMeans(M)**2)*nrow(M)/(nrow(M)-1)
    return(result)
}


# may be useful when applying multiple contrasts, but we need
# a solution to apply contrasts and deduce ranks rapidly
order_sparse = function(M) {
    R = tibble(values = M@x,
               row_index = M@i,
               column_index = rep(1:ncol(M), diff(M@p))) %>%
        arrange(column_index, values)
    return(R)
}


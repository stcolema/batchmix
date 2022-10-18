#' @title Generate batch data
#' @description Generate data from K multivaraite normal or multivariate t
#' distributions with additional noise from batches. Assumes independence across
#' columns. In each column the parameters are randomly permuted for both the
#' groups and batches.
#' @param N The number of items (rows) to generate.
#' @param P The number of columns in the generated dataset.
#' @param group_rates A vector of the group rates for the classes within a column.
#' @param batch_rates A vector of the batch rates for the classes within a column.
#' This is used to create a variable which has the sum of the appropriate batch 
#' and class rate, it might be better interpreted as the batch effect on the 
#' observed rate.
#' @param group_weights One of either a K x B matrix of the expected proportion 
#' of each batch in each group or a K-vector of the expected proportion of the
#' entire dataset in each group.
#' @param batch_weights A vector of the expected proportion of N in each batch.
#' @param frac_known The number of items with known labels.
#' @param permute_variables Logical indicating if group and batch means and
#' standard deviations should be permuted in each column or not (defaults to
#' ``TRUE``).
#' @param scale_data Logical indicating if data should be mean centred and 
#' standardised (defaults to ``FALSE``).
#' @return A list of 5 objects; the data generated from the groups with and
#' without batch effects, the label indicating the generating group, the
#' batch label and the vector indicating training versus test.
#' @export
generateBatchDataLogPoisson <- function(N,
                                        P,
                                        group_rates,
                                        batch_rates,
                                        group_weights,
                                        batch_weights,
                                        frac_known = 0.2,
                                        permute_variables = TRUE,
                                        scale_data = FALSE
    ) {
  
  
  # The number of batches and groups to generate
  B <- length(batch_rates)
  K <- length(group_rates)
  
  # Allow for varying groups within batches
  varying_group_within_batch <- is.matrix(group_weights)
  
  # The batch labels for the N points
  batch_IDs <- sample(seq(1, B), N, replace = TRUE, prob = batch_weights)
  
  # Generate group membership, potentially allowing imbalance across batches
  group_IDs <- generateGroupIDsInSimulator(N,
    K,
    B,
    batch_IDs,
    group_weights,
    varying_group_within_batch
  )
  
  # Fixed labels
  fixed <- sample(seq(0, 1), N,
    replace = TRUE,
    prob = c(1 - frac_known, frac_known)
  )
  
  classes <- unique(group_IDs)
  batches <- unique(batch_IDs)
  observed_data <- true_data <- matrix(nrow = N, ncol = P)
  
  lambdas_pk <- group_rates
  lambdas_pb <- batch_rates
  
  for(p in seq(1, P)) {
    if(permute_variables) {
      lambdas_pk <- sample(group_rates)
      lambdas_pb <- sample(batch_rates)
    }
    for(b in batches) {
      batch_inds <- which(batch_IDs == b)
      for(k in classes) {
        class_in_batch_inds <- which((group_IDs == k) & (batch_IDs == b))
        n_kb <- length(class_in_batch_inds)
        if(n_kb > 0) {
          batch_free_count <- rpois(n_kb, lambdas_pk[k])
          batch_effect_on_count <- rpois(n_kb, lambdas_pb[b])
          observed_count <- batch_free_count + batch_effect_on_count
          normal_noise <- rnorm(n_kb)
          observed_data[class_in_batch_inds, p] <- log(1.0 + observed_count) + normal_noise
          true_data[class_in_batch_inds, p] <- log(1.0 + batch_free_count) + normal_noise
        }
      }
    }
  }
  row.names(true_data) <- row.names(observed_data) <- paste0("Person_", seq(1, N))
  
  if(scale_data) {
    observed_data <- scale(observed_data)
    true_data <- scale(true_data)
  }
  
  # Return the data, the data without batch effects, the allocation labels and
  # the batch labels.
  list(
    observed_data = observed_data,
    corrected_data = true_data,
    group_IDs = group_IDs,
    batch_IDs = batch_IDs,
    fixed = fixed
  )
}


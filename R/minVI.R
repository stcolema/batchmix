#' @title Minimium VI
#' @description Local implementation of S. Wade's `minVI` function from their
#' `mcclust.ext` package (available from github).
#' Reimplemented here to avoid dependency on a non-CRAN package and we have
#' dropped the `greedy` method. Finds the optimal partition by minimising the
#' lower bound to the Variation of Information obtained from Jensen's inequality
#' where the expectation and log are reversed. For full details please see the
#' aforementioned package and Wade and Ghahramani, 2018, 'Bayesian Cluster
#' Analysis: Point Estimation and Credible Balls (with Discussion)'.
#' @param psm The posterior similarity matrix for a set of clustering MCMC
#' samples such as is returned by the `createSimilarityMat` function.
#' @param cls.draw The set of clustering MCMC samples used to generate `psm`.
#' Only required if `method` is one of `'draws'` or `'all'`.
#' @param method String indicating which method is used to find the point
#' estimate clustering. Must be one of `'avg'`, `'comp'`, `'draws'` or `'all'`.
#' Defaults to `'avg'`. If `'all'` is passed the three methods are all applied
#' to return different choices of point clustering.
#' @param max.k The maximum number of clusters to consider. Only used by the
#' `'comp'` and `'avg'` methods. Defaults to one-quarter the number of data
#' points rounded up.
#' @returns If `method` is `'all'` returns a matrix of four clusterings, one for
#' each method and a repeat of that which performs best based on minimising the
#' Variation of Information between the clustering and the PSM. Otherwise
#' returns a vector. This is annotated with the attribute `"info"`, a named list
#' describing:
#'
#' * `.$loss`: the loss score used (Variation of Information)
#'
#' * `.$maxNClusters`: the `max.k` value used by the `'comp'` and `'avg'` methods
#'
#' * `.$expectedLoss`: the estimated minimum Variation of Information for the point
#' clustering(s)
#'
#' * `.$method`: the point method used to infer the clustering(s)
#'
#' Names are due to legacy reasons - this function is replacing the
#' `salso::salso` function and name choices are to minimise workflow damage.
#' @examples
#' \dontrun{
#' # MCMC samples and BIC vector
#' mcmc_outputs <- runMCMCChains(
#'   X,
#'   n_chains,
#'   R,
#'   thin,
#'   batch_vec,
#'   type
#' )
#'
#' # Note that in this toy example we have not applied a burn in
#' psm <- createSimilarityMat(mcmc_outputs[[1]]$samples)
#' cl_est <- minVI(psm, mcmc_outputs[[1]]$samples)
#' }
#' @importFrom stats hclust as.dist cutree
#' @export
minVI <- function(psm,
                  cls.draw = NULL,
                  method = "avg",
                  max.k = NULL) {
  wrong_method <- !(method %in% c("avg", "comp", "draws", "all"))
  if (wrong_method) {
    stop("method must be one of 'avg', 'comp', 'draws' or 'all'.")
  }

  method <- match.arg(method, choices = method)
  if (method %in% c("draws", "all") & is.null(cls.draw)) {
    stop("cls.draw must be provided if method=''draws''")
  }

  # If no maximum number of clusters is passed, set this to quarter of the
  # number of datapoints - note this is only used by 'avg' and 'comp' methods
  if (is.null(max.k)) {
    max.k <- ceiling(dim(psm)[1] / 4)
    k_inds <- seq(1, max.k)
  }

  # If using the 'average' method
  if (method == "avg" | method == "all") {
    hclust.avg <- hclust(as.dist(1 - psm), method = "average")
    cls.avg <- t(apply(matrix(k_inds), 1, function(x) cutree(hclust.avg, k = x)))
    VI.avg <- VI.lb(cls.avg, psm)
    val.avg <- min(VI.avg)
    cl.avg <- cls.avg[which.min(VI.avg), ]
    if (method == "avg") {
      avg_estimate <- cl.avg
      avg_info <- list(
        loss = "VI",
        maxNClusters = max.k,
        expectedLoss = val.avg,
        method = "avg"
      )
      attr(avg_estimate, "info") <- avg_info
      return(avg_estimate)
    }
  }

  if (method == "comp" | method == "all") {
    # if (is.null(max.k)) max.k <- ceiling(dim(psm)[1] / 8)
    hclust.comp <- hclust(as.dist(1 - psm), method = "complete")
    cls.comp <- t(apply(matrix(k_inds), 1, function(x) cutree(hclust.comp, k = x)))
    VI.comp <- VI.lb(cls.comp, psm)
    val.comp <- min(VI.comp)
    cl.comp <- cls.comp[which.min(VI.comp), ]
    if (method == "comp") {
      comp_estimate <- cl.comp
      comp_info <- list(
        loss = "VI",
        maxNClusters = max.k,
        expectedLoss = val.comp,
        method = "comp"
      )
      attr(comp_estimate, "info") <- comp_info
      return(comp_estimate)
    }
  }

  if (method == "draws" | method == "all") {
    n <- ncol(psm)
    EVI_lb_local <- function(c) {
      f <- 0
      n_inds <- seq(1, n)
      for (i in n_inds) {
        ind <- (c == c[i])
        f <- f + (log2(sum(ind))
        + log2(sum(psm[i, ]))
          - 2 * log2(sum(ind * psm[i, ]))
        ) / n
      }
      return(f)
    }
    VI.draws <- apply(cls.draw, 1, EVI_lb_local)
    val.draws <- min(VI.draws)
    cl.draw <- cls.draw[which.min(VI.draws), ]
    names(cl.draw) <- NULL
    if (method == "draws") {
      draws_estimate <- cl.draw
      draws_info <- list(
        loss = "VI",
        maxNClusters = max.k,
        expectedLoss = val.draws,
        method = "draws"
      )
      attr(draws_estimate, "info") <- draws_info
      return(draws_estimate)
    }
  }

  # If using "all" compile the information
  vals <- c(val.avg, val.comp, val.draws)
  cls <- rbind(cl.avg, cl.comp, cl.draw)
  cls <- rbind(cls[which.min(vals), ], cls)
  vals <- c(min(vals), vals)
  rownames(cls) <- names(vals) <- c("best", "avg", "comp", "draws")
  colnames(cls) <- NULL

  estimate <- cls
  info <- list(
    loss = "VI",
    maxNClusters = max.k,
    expectedLoss = vals,
    method = "all"
  )
  attr(estimate, "info") <- info
  estimate
}

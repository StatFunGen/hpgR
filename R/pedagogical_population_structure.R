#' Population structure and ancestry likelihoods
#'
#' PCA centers on 2p and scales by sqrt(2p(1-p)); markers with expected genotype standard deviation at most 0.05 are excluded. The admixture EM uses independent binomial genotype likelihoods and k-means initialization on up to four supplied PC scores.
#'
#' @name models_structure
#' @param genotype Complete numeric genotype counts 0,1,2: a people-by-marker matrix for PCA and admixture, or a vector for one person in ancestryLikelihood.
#' @param components Positive integer maximum number of PCs to retain, limited by the data rank.
#' @param k Positive integer number of ancestry sources, at most the number of people and distinct score rows.
#' @param scores Finite numeric PC-score matrix with the same people in rows as genotype.
#' @param max_iter Positive integer EM iteration limit.
#' @param per_cell_tol Positive per-genotype-cell log likelihood change tolerance.
#' @param q_tol Positive maximum ancestry-proportion change tolerance; both tolerances must hold for ten consecutive updates.
#' @param q Source probability vector for one person, summing to one.
#' @param frequency Marker-by-source matrix of reference allele frequencies in [0,1].
#' @param log Logical; return the log likelihood when TRUE.
#' @returns genotypePCA returns scores, allele_frequency, keep, variance_fraction and rotation (retained markers by PCs). fitAdmixture returns q (people by sources), p (markers by sources), log_likelihood, converged and iterations, warning on an exhausted limit. ancestryLikelihood returns a likelihood or log likelihood for one person.
#' @examples
#' set.seed(42)
#' G <- simulateGenotypes(20, rep(.3, 10))
#' pc <- genotypePCA(G, components=2)
#' dim(pc$scores)
#' ancestryLikelihood(c(.4,.6), c(0,1), cbind(c(.2,.8),c(.7,.3)))
NULL

#' @rdname models_structure
#' @export
ancestryLikelihood <- function(q, genotype, frequency, log = FALSE) {
    checkInputs(genotype = dosageVector,
                frequency = list(numericMatrix, minimum = 0, maximum = 1),
                q = list(probability, length = ncol(frequency)), log = logicalScalar)
    if (nrow(frequency) != length(genotype) || abs(sum(q) - 1) >= 1e-10)
        stop("frequency must have one row per marker and q must sum to 1", call. = FALSE)
    f <- as.vector(frequency %*% q)
    if (log) sum(dbinom(genotype, size = 2, prob = f, log = TRUE)) else
        prod(dbinom(genotype, size = 2, prob = f))
}

#' @rdname models_structure
#' @export
genotypePCA <- function(genotype, components = 4L) {
  checkInputs(genotype = list(dosageMatrix, rows = 2L),
              components = list(scalarInteger, minimum = 1L))
  allele_frequency <- colMeans(genotype) / 2
  marker_sd <- sqrt(2 * allele_frequency * (1 - allele_frequency))
  ## Exclude nearly monomorphic SNPs with an unstable scaling denominator.
  keep <- marker_sd > 0.05
  if (!any(keep)) stop("genotype has no markers with sufficient variation for PCA", call. = FALSE)
  standardized <- standardizeGenotypes(genotype[, keep, drop = FALSE],
                                      allele_frequency[keep])
  fit <- prcomp(standardized, center = FALSE, scale. = FALSE,
                rank. = components)
  variance_fraction <- fit$sdev^2 / sum(fit$sdev^2)
  list(scores = fit$x, allele_frequency = allele_frequency, keep = keep,
       variance_fraction = variance_fraction[seq_len(ncol(fit$x))], rotation = fit$rotation)
}

#' Project genotypes onto fitted reference principal components
#'
#' Apply the reference allele frequencies, retained markers and loadings from
#' genotypePCA. Genotypes have the same marker order as the reference data.
#' @param genotype Complete people-by-marker genotype matrix, counting allele 1.
#' @param pca Fitted list returned by genotypePCA.
#' @return A people-by-PC matrix of projected scores.
#' @examples
#' G <- rbind(c(0,0,1), c(1,1,2), c(2,2,0), c(0,1,1))
#' pc <- genotypePCA(G, 2)
#' projectGenotypePCA(G[1,,drop=FALSE], pc)
#' @export
projectGenotypePCA <- function(genotype, pca) {
    checkInputs(genotype = dosageMatrix)
    if (!is.list(pca) || !is.logical(pca$keep) || anyNA(pca$keep) ||
        length(pca$keep) != ncol(genotype) || !any(pca$keep))
        stop("pca must retain markers in the supplied genotype order", call. = FALSE)
    probability(pca$allele_frequency, "pca$allele_frequency", ncol(genotype))
    numericMatrix(pca$rotation, "pca$rotation")
    if (nrow(pca$rotation) != sum(pca$keep))
        stop("pca loadings must match the retained markers", call. = FALSE)
    if (!is.null(colnames(genotype)) && !is.null(rownames(pca$rotation)) &&
        !identical(colnames(genotype)[pca$keep], rownames(pca$rotation)))
        stop("genotype marker names must match the fitted marker order", call. = FALSE)
    standardizeGenotypes(genotype[, pca$keep, drop = FALSE],
                         pca$allele_frequency[pca$keep]) %*% pca$rotation
}

#' Estimate ancestry fractions using two fixed source panels
#'
#' For each person, maximize ancestryLikelihood over the fraction q from the
#' second source, with marker probability (1-q)*frequency[,1]+q*frequency[,2].
#' The likelihood treats marker genotypes as independent binomial counts.
#' @param genotype Complete people-by-marker genotype matrix, counting allele 1.
#' @param frequency Marker-by-two matrix of supplied source allele frequencies.
#' @param tolerance Positive tolerance for stats::optimize.
#' @return A vector of fitted second-source fractions, one per person.
#' @examples
#' F <- cbind(c(.2,.3,.1), c(.8,.7,.9))
#' fitSourceAncestry(rbind(c(0,1,0), c(2,1,2)), F)
#' @export
fitSourceAncestry <- function(genotype, frequency, tolerance = 1e-6) {
    checkInputs(genotype = dosageMatrix,
                frequency = list(numericMatrix, minimum = 0, maximum = 1),
                tolerance = positiveScalar)
    if (ncol(frequency) != 2L || nrow(frequency) != ncol(genotype))
        stop("frequency must give two source frequencies per marker", call. = FALSE)
    if (all(frequency[, 1L] == frequency[, 2L]))
        stop("source frequencies must differ to identify ancestry", call. = FALSE)
    apply(genotype, 1L, function(g) {
        objective <- function(q) -ancestryLikelihood(c(1 - q, q), g, frequency, log = TRUE)
        if (!is.finite(objective(.5)))
            stop("genotypes have zero probability under the source panels", call. = FALSE)
        fit <- stats::optimize(objective, c(0, 1), tol = tolerance)
        candidates <- c(0, fit$minimum, 1)
        candidates[which.min(vapply(candidates, objective, numeric(1L)))]
    })
}

#' @rdname models_structure
#' @export
fitAdmixture <- function(genotype, k, scores, max_iter = 2000L,
                         per_cell_tol = 1e-7, q_tol = 1e-4) {
  checkInputs(genotype = list(dosageMatrix, rows = 2L),
              k = list(scalarInteger, minimum = 1L, maximum = nrow(genotype)),
              scores = list(numericMatrix, rows = 2L),
              max_iter = list(scalarInteger, minimum = 1L), per_cell_tol = positiveScalar,
              q_tol = positiveScalar)
  if (nrow(scores) != nrow(genotype))
    stop("scores must have one row per person", call. = FALSE)
  n <- nrow(genotype)
  m <- ncol(genotype)
  centers <- kmeans(scores[, seq_len(min(ncol(scores), 4L)), drop = FALSE],
                    centers = k, nstart = 10L)$cluster
  q <- matrix(0.1 / k, nrow = n, ncol = k)
  q[cbind(seq_len(n), centers)] <- 0.9 + 0.1 / k
  p <- vapply(seq_len(k), function(j) {
    colSums(sweep(genotype, 1, q[, j], "*")) / (2 * sum(q[, j]))
  }, numeric(m))
  p <- matrix(p, nrow = m, ncol = k)
  p[p < 0.005] <- 0.005
  p[p > 0.995] <- 0.995

  log_likelihood <- numeric(max_iter)
  stable_updates <- 0L
  for (iteration in seq_len(max_iter)) {
    previous_q <- q
    f <- q %*% t(p)
    f[f < 0.001] <- 0.001
    f[f > 0.999] <- 0.999
    copy_count <- matrix(0, n, k)
    allele_count <- matrix(0, m, k)
    source_count <- matrix(0, m, k)

    for (j in seq_len(k)) {
      alt <- sweep(matrix(p[, j], n, m, byrow = TRUE), 1, q[, j], "*") / f
      ref <- sweep(matrix(1 - p[, j], n, m, byrow = TRUE), 1,
                   q[, j], "*") / (1 - f)
      assigned_alt <- genotype * alt
      assigned_ref <- (2 - genotype) * ref
      assigned <- assigned_alt + assigned_ref
      copy_count[, j] <- rowSums(assigned)
      allele_count[, j] <- colSums(assigned_alt)
      source_count[, j] <- colSums(assigned)
    }
    ## Keep a nearly unused source from disappearing numerically.
    q <- pmax(copy_count / (2 * m), 1e-6)
    q <- q / rowSums(q)
    p <- allele_count / source_count
    p[p < 0.005] <- 0.005
    p[p > 0.995] <- 0.995
    ## Record the likelihood after the update, using the returned P and Q.
    log_likelihood[iteration] <- sum(dbinom(genotype, 2, q %*% t(p), log = TRUE))
    if (iteration > 1L) {
      gain <- log_likelihood[iteration] - log_likelihood[iteration - 1L]
      small_change <- abs(gain) / (n * m) < per_cell_tol &&
                      max(abs(q - previous_q)) < q_tol
      stable_updates <- if (small_change) stable_updates + 1L else 0L
      if (stable_updates >= 10L) break
    }
  }
  stopifnot(max(abs(rowSums(q) - 1)) < 1e-8,
            all(p >= 0), all(p <= 1))
  converged <- stable_updates >= 10L
  if (!converged) warning("Iteration limit reached; inspect the trace or increase max_iter.")
  list(q = q, p = p, log_likelihood = log_likelihood[seq_len(iteration)],
       converged = converged, iterations = iteration)
}

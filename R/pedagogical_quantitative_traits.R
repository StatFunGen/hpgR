#' Additive quantitative-trait and genotype covariance models
#'
#' Genotype simulation uses independent binomial marker counts. Standardization uses the supplied marker frequencies. The polygenic simulator rescales realized genetic values to sample variance V_G. The profile likelihood assumes a centered phenotype and covariance proportional to h2 K + (1-h2) I; it profiles total variance by maximum likelihood.
#'
#' @name models_traits
#' @param label One character population label.
#' @param mean_height Finite population phenotype mean.
#' @param additive_variance Nonnegative additive genetic variance.
#' @param environmental_variance Nonnegative environmental variance.
#' @param n Positive integer number of people.
#' @param standardized_genotype Finite numeric matrix of standardized genotypes, with people in rows and markers in columns.
#' @param V_G Nonnegative target sample genetic variance; at least two people are required.
#' @param V_E Nonnegative independent environmental variance.
#' @param n_pairs Positive integer number of sibling pairs.
#' @param V_A Nonnegative additive genetic variance.
#' @param V_C Nonnegative shared-environment variance.
#' @param genotype Complete numeric people-by-marker genotype-count matrix (0,1,2).
#' @param frequency One allele frequency per marker; standardization requires values strictly between zero and one.
#' @param h2 Scalar heritability in [0,1].
#' @param eigenvalues Nonnegative eigenvalues of the genetic relationship matrix.
#' @param rotated_phenotype Centered phenotype rotated into the corresponding eigenvectors; same length as eigenvalues.
#' @returns simulateGenotypes returns a people-by-marker matrix, including singleton dimensions. standardizeGenotypes returns a matrix of the same shape. geneticRelationship returns a people-by-people matrix. geneticProfileNLL returns a scalar objective (Inf for singular covariance). simulatePolygenicTrait returns marker_effect, genetic_value, environmental_deviation and centered phenotype, plus drawn_effect and effect_scale; marker_effect is the rescaled effect satisfying Z times marker_effect = genetic_value. simulateSiblingTraits returns realized sharing and both phenotypes. simulateGaussianTrait returns a data frame of population, genetic_value, environmental_deviation and height.
#' @examples
#' set.seed(42)
#' G <- simulateGenotypes(20, rep(.3,5))
#' Z <- standardizeGenotypes(G, rep(.3,5))
#' trait <- simulatePolygenicTrait(Z, .5, .5)
#' stopifnot(isTRUE(all.equal(drop(Z %*% trait$marker_effect), trait$genetic_value)))
NULL

#' @rdname models_traits
#' @export
simulateGenotypes <- function(n, frequency) {
    .scalarInteger(n, "n", 1L)
    .probability(frequency, "frequency", length = NULL)
    result <- vapply(frequency, function(p) rbinom(n, size = 2, prob = p), numeric(n))
    matrix(result, nrow = n, dimnames = list(NULL, names(frequency)))
}

# Continuous additive trait conditional on realized sibling sharing.
#' @rdname models_traits
#' @export
simulateSiblingTraits <- function(n_pairs, V_A, V_C, V_E, mean_height = 0) {
    .scalarInteger(n_pairs, "n_pairs", 1L)
    .positiveScalar(V_A, "V_A", zero = TRUE)
    .positiveScalar(V_C, "V_C", zero = TRUE)
    .positiveScalar(V_E, "V_E", zero = TRUE)
    .numericVector(mean_height, "mean_height", 1L)
    ibd_sharing <- rbeta(n_pairs, shape1 = 50, shape2 = 50)
    z1 <- rnorm(n_pairs)
    z2 <- rnorm(n_pairs)
    genetic_1 <- sqrt(V_A) * z1
    genetic_2 <- sqrt(V_A) * (ibd_sharing * z1 + sqrt(1 - ibd_sharing^2) * z2)
    shared_environment <- rnorm(n_pairs, 0, sqrt(V_C))
    environment_1 <- rnorm(n_pairs, 0, sqrt(V_E))
    environment_2 <- rnorm(n_pairs, 0, sqrt(V_E))
    phenotype_1 <- mean_height + genetic_1 + shared_environment + environment_1
    phenotype_2 <- mean_height + genetic_2 + shared_environment + environment_2
    list(ibd_sharing = ibd_sharing, phenotype_1 = phenotype_1, phenotype_2 = phenotype_2)
}

#' @rdname models_traits
#' @export
standardizeGenotypes <- function(genotype, frequency = colMeans(genotype) / 2) {
    .dosageMatrix(genotype)
    .probability(frequency, "frequency", ncol(genotype))
    if (any(frequency == 0 | frequency == 1))
        stop("frequency must be strictly between 0 and 1 for standardization", call. = FALSE)
    centered <- sweep(genotype, 2, 2 * frequency, '-')
    sweep(centered, 2, sqrt(2 * frequency * (1 - frequency)), '/')
}

#' @rdname models_traits
#' @export
geneticRelationship <- function(standardized_genotype) {
    .numericMatrix(standardized_genotype, "standardized_genotype")
    tcrossprod(standardized_genotype) / ncol(standardized_genotype)
}

# Profile ML with a centered phenotype and one genetic covariance component.
# This is the classroom likelihood; it is not a REML implementation.
#' @rdname models_traits
#' @export
geneticProfileNLL <- function(h2, eigenvalues, rotated_phenotype) {
    .probability(h2, "h2")
    .numericVector(eigenvalues, "eigenvalues", minimum = 0)
    .numericVector(rotated_phenotype, "rotated_phenotype", length(eigenvalues))
    if (!any(rotated_phenotype != 0))
        stop("rotated_phenotype must have positive variation", call. = FALSE)
    n_people <- length(rotated_phenotype)
    covariance_eigenvalue <- h2 * eigenvalues + (1 - h2)
    if (any(covariance_eigenvalue <= 0)) return(Inf)
    total_variance <- mean(rotated_phenotype^2 / covariance_eigenvalue)
    0.5 * (sum(log(covariance_eigenvalue)) + n_people * log(total_variance) +
        sum(rotated_phenotype^2 / (total_variance * covariance_eigenvalue)))
}

#' @rdname models_traits
#' @export
simulatePolygenicTrait <- function(standardized_genotype, V_G, V_E) {
    .numericMatrix(standardized_genotype, "standardized_genotype", rows = 2L)
    .positiveScalar(V_G, "V_G", zero = TRUE)
    .positiveScalar(V_E, "V_E", zero = TRUE)
    n_snps <- ncol(standardized_genotype)
    n_people <- nrow(standardized_genotype)
    marker_effect <- rnorm(n_snps, 0, sqrt(V_G / n_snps))
    genetic_value <- as.vector(standardized_genotype %*% marker_effect)
    drawn_effect <- marker_effect
    if (V_G > 0 && var(genetic_value) == 0)
        stop("standardized_genotype must support positive genetic variance", call. = FALSE)
    effect_scale <- if (V_G == 0) 0 else sqrt(V_G / var(genetic_value))
    genetic_value <- genetic_value * effect_scale
    marker_effect <- marker_effect * effect_scale
    environmental_deviation <- rnorm(n_people, 0, sqrt(V_E))
    phenotype <- genetic_value + environmental_deviation
    phenotype <- phenotype - mean(phenotype)
    list(marker_effect = marker_effect, genetic_value = genetic_value,
        environmental_deviation = environmental_deviation, phenotype = phenotype,
        drawn_effect = drawn_effect, effect_scale = effect_scale)
}

#' @rdname models_traits
#' @export
simulateGaussianTrait <- function(label, mean_height, additive_variance,
                            environmental_variance, n) {
  .scalarInteger(n, "n", 1L)
  .numericVector(mean_height, "mean_height", 1L)
  .positiveScalar(additive_variance, "additive_variance", zero = TRUE)
  .positiveScalar(environmental_variance, "environmental_variance", zero = TRUE)
  if (!is.character(label) || length(label) != 1L || is.na(label))
    stop("label must be one character string", call. = FALSE)
  additive_value <- rnorm(n, 0, sqrt(additive_variance))
  environmental_deviation <- rnorm(n, 0, sqrt(environmental_variance))

  data.frame(
    population = label,
    genetic_value = additive_value,
    environmental_deviation = environmental_deviation,
    height = mean_height + additive_value + environmental_deviation
  )
}

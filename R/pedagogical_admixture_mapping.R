# Local-ancestry association conditional on global ancestry.

#' @rdname models_structure
#' @export
simulateAncestryAssociation <- function(people, loci, target, intercept = -2.5,
                                       local_effect = 1.4, global_effect = 0.4) {
    scalarInteger(people, "people", 1L)
    scalarInteger(loci, "loci", 1L)
    scalarInteger(target, "target", 1L, loci)
    numericVector(c(intercept, local_effect, global_effect), "regression coefficients", 3L)
    global <- runif(people, 0.2, 0.8)
    local <- matrix(vapply(seq_len(loci), function(j) rbinom(people, 2, global),
                           numeric(people)), nrow = people)
    disease_probability <- plogis(intercept + local_effect * local[, target] +
        global_effect * global)
    disease <- rbinom(people, 1, disease_probability)
    list(global = global, local = local, disease_probability = disease_probability,
        disease = disease)
}

#' @rdname models_structure
#' @export
fitAncestryAssociation <- function(disease, local, global) {
    dosageMatrix(local, "local")
    probability(disease, "disease", nrow(local))
    if (any(!disease %in% 0:1)) stop("disease must contain 0 or 1", call. = FALSE)
    probability(global, "global", nrow(local))
    t(vapply(seq_len(ncol(local)), function(j) {
        data <- data.frame(disease = disease, local = local[, j], global = global)
        fit <- glm(disease ~ local + global, data = data, family = binomial())
        c(estimate = coef(fit)['local'], se = summary(fit)$coefficients['local', 'Std. Error'])
    }, numeric(2)))
}

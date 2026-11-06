# Shared input checks. A list supplies additional arguments to a validator.
# checkInputs(n = scalarInteger, p = list(probability, length = NULL))
checkInputs <- function(...) {
    rules <- as.list(substitute(list(...)))[-1L]
    labels <- names(rules)
    if (length(rules) && (is.null(labels) || any(!nzchar(labels)) ||
        anyDuplicated(labels)))
        stop("checkInputs requires unique argument names", call. = FALSE)
    caller <- parent.frame()
    for (label in labels) {
        rule <- eval(rules[[label]], envir = caller)
        if (is.function(rule)) rule <- list(rule)
        if (!is.list(rule) || !length(rule) || !is.function(rule[[1L]]))
            stop("invalid input check for ", label, call. = FALSE)
        arguments <- c(list(get(label, envir = caller, inherits = FALSE)), rule[-1L])
        if ("name" %in% names(formals(rule[[1L]]))) arguments$name <- label
        do.call(rule[[1L]], arguments)
    }
    invisible(NULL)
}

scalarInteger <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x) || x > maximum || x > .Machine$integer.max) {
        stop(name, " must be an integer >= ", minimum, call. = FALSE)
    }
    as.integer(x)
}

positiveScalar <- function(x, name, zero = FALSE) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
        if (zero) x < 0 else x <= 0) {
        stop(name, " must be a finite ", if (zero) "nonnegative" else "positive",
             " number", call. = FALSE)
    }
    x
}

numericVector <- function(x, name, length = NULL, minimum = -Inf, maximum = Inf,
                           empty = FALSE, integer = FALSE, nondecreasing = FALSE) {
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
        (!empty && !base::length(x)) || any(!is.finite(x)) ||
        any(x < minimum | x > maximum) ||
        (!is.null(length) && base::length(x) != length)) {
        stop(name, " must be a finite numeric vector",
             if (!is.null(length)) paste0(" of length ", length) else "",
             if (is.finite(minimum) || is.finite(maximum))
                 paste0(" with values in [", minimum, ", ", maximum, "]") else "",
             call. = FALSE)
    }
    if (integer && any(x != floor(x)))
        stop(name, " must contain integers", call. = FALSE)
    if (nondecreasing && any(diff(x) < 0))
        stop(name, " must be nondecreasing", call. = FALSE)
    x
}

probability <- function(x, name, length = 1L, open = FALSE, sumToOne = FALSE) {
    numericVector(x, name, length, 0, 1)
    if (open && any(x == 0 | x == 1))
        stop(name, " must be strictly between 0 and 1", call. = FALSE)
    if (sumToOne && abs(sum(x) - 1) >= 1e-10)
        stop(name, " must sum to 1", call. = FALSE)
    x
}

numericMatrix <- function(x, name, rows = 1L, columns = 1L,
                          minimum = -Inf, maximum = Inf) {
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        nrow(x) < rows || ncol(x) < columns || any(!is.finite(x))) {
        stop(name, " must be a finite numeric matrix with at least ", rows,
             " row(s) and ", columns, " column(s)", call. = FALSE)
    }
    if (any(x < minimum | x > maximum))
        stop(name, " must have values in [", minimum, ", ", maximum, "]",
             call. = FALSE)
    x
}

binaryVector <- function(x, name, length = NULL, missing = FALSE, logical = FALSE) {
    if (!missing) probability(x, name, length) else {
        if (!(is.numeric(x) || (logical && is.logical(x))) || is.complex(x) ||
            !is.null(dim(x)) || (!is.null(length) && base::length(x) != length))
            stop(name, " must give one binary allele or NA per marker", call. = FALSE)
    }
    if (any(!is.na(x) & !x %in% 0:1))
        stop(name, " must contain binary alleles 0 or 1", call. = FALSE)
    x
}

binaryMatrix <- function(x, name) {
    numericMatrix(x, name)
    if (any(!x %in% 0:1)) stop(name, " must contain 0 or 1", call. = FALSE)
    x
}

logicalScalar <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
        stop(name, " must be TRUE or FALSE", call. = FALSE)
    x
}

transitionProbabilities <- function(x, name, states, observations) {
    varying <- length(dim(x)) == 3L
    expected <- if (varying) c(states, states, observations - 1L) else c(states, states)
    if (!is.numeric(x) || is.complex(x) ||
        !identical(as.integer(dim(x)), as.integer(expected)) ||
        any(!is.finite(x)) || any(x < 0))
        stop(name, " dimensions or probabilities are invalid", call. = FALSE)
    if (varying && observations > 1L) {
        for (j in seq_len(observations - 1L))
            if (any(abs(rowSums(matrix(x[, , j], states, states)) - 1) > 1e-10))
                stop(name, " rows must sum to one", call. = FALSE)
    } else if (!varying && any(abs(rowSums(x) - 1) > 1e-10)) {
        stop(name, " rows must sum to one", call. = FALSE)
    }
    x
}

dosageVector <- function(x, name) {
    numericVector(x, name, minimum = 0, maximum = 2)
    if (any(x != floor(x))) stop(name, " must contain genotype counts 0, 1 or 2",
                                call. = FALSE)
    x
}

dosageMatrix <- function(x, name = "genotype", rows = 1L) {
    numericMatrix(x, name, rows)
    if (any(!x %in% 0:2)) stop(name, " must contain genotype counts 0, 1 or 2",
                             call. = FALSE)
    x
}

sameLength <- function(x, y, names) {
    if (length(x) != length(y)) stop(paste(names, collapse = " and "),
        " must have the same length", call. = FALSE)
    invisible(NULL)
}

fitnessWeights <- function(w) {
    numericVector(w, "w", length = 3L, minimum = 0)
    if (!any(w > 0)) stop("w must include positive fitness", call. = FALSE)
    w
}

# A cleanup closure also restores sessions that had no RNG state on entry.
rngRestorer <- function() {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    old_kind <- RNGkind()
    function() {
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
        invisible(NULL)
    }
}

genotypeMatrix <- function(x) {
    if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) ||
        any(!is.finite(x)) || any(x < 0 | x != floor(x)) ||
        any(x > .Machine$integer.max)) {
        stop("genotypes must be a nonempty matrix of nonnegative integer alleles",
             call. = FALSE)
    }
    storage.mode(x) <- "integer"
    x
}

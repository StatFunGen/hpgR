# Shared input checks for the course's R interfaces and classroom examples.
.scalarInteger <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x) || x > maximum || x > .Machine$integer.max) {
        stop(name, " must be an integer >= ", minimum, call. = FALSE)
    }
    as.integer(x)
}

.positiveScalar <- function(x, name, zero = FALSE) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
        if (zero) x < 0 else x <= 0) {
        stop(name, " must be a finite ", if (zero) "nonnegative" else "positive",
             " number", call. = FALSE)
    }
    x
}

.numericVector <- function(x, name, length = NULL, minimum = -Inf, maximum = Inf,
                           empty = FALSE) {
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
    x
}

.probability <- function(x, name, length = 1L) {
    .numericVector(x, name, length, 0, 1)
}

.numericMatrix <- function(x, name, rows = 1L, columns = 1L) {
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        nrow(x) < rows || ncol(x) < columns || any(!is.finite(x))) {
        stop(name, " must be a finite numeric matrix with at least ", rows,
             " row(s) and ", columns, " column(s)", call. = FALSE)
    }
    x
}

.dosageVector <- function(x, name) {
    .numericVector(x, name, minimum = 0, maximum = 2)
    if (any(x != floor(x))) stop(name, " must contain genotype counts 0, 1 or 2",
                                call. = FALSE)
    x
}

.dosageMatrix <- function(x, name = "genotype", rows = 1L) {
    .numericMatrix(x, name, rows)
    if (any(!x %in% 0:2)) stop(name, " must contain genotype counts 0, 1 or 2",
                             call. = FALSE)
    x
}

.sameLength <- function(x, y, names) {
    if (length(x) != length(y)) stop(paste(names, collapse = " and "),
        " must have the same length", call. = FALSE)
    invisible(NULL)
}

.fitnessWeights <- function(w) {
    .numericVector(w, "w", length = 3L, minimum = 0)
    if (!any(w > 0)) stop("w must include positive fitness", call. = FALSE)
    w
}

# A cleanup closure also restores sessions that had no RNG state on entry.
.rngRestorer <- function() {
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

.genotypeMatrix <- function(x) {
    if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) ||
        any(!is.finite(x)) || any(x < 0 | x != floor(x)) ||
        any(x > .Machine$integer.max)) {
        stop("genotypes must be a nonempty matrix of nonnegative integer alleles",
             call. = FALSE)
    }
    storage.mode(x) <- "integer"
    x
}

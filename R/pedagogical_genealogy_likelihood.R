#' Sequence constraints and mutation-clock fits for four copies
#'
#' An error-free nonrecombining region follows one rooted genealogy with infinite-sites mutation. Shared derived variants constrain clades. Compatible ranked binary trees are fitted by a conditional Poisson likelihood; equal best fits retain unresolved ages. Linked constraints report minimal compatible blocks and possible cut gaps.
#'
#' @name models_tree_inference
#' @param G Numeric binary haplotype matrix, sampled copies in rows and polarized sites in columns; zero is ancestral. Most pti_ functions require exactly four rows.
#' @param mu Positive mutation rate per base per generation.
#' @param L Positive analyzed sequence length in bases.
#' @param tie_tolerance Nonnegative log likelihood tolerance for retaining tied best fits.
#' @param tree Recursively nested list of labelled tips.
#' @param theta Nonnegative interval mutation exposures (mutation rate times length times generation durations).
#' @param A Numeric edge-by-ranked-interval design matrix from pti_rank_designs.
#' @param counts Nonnegative derived-variant counts assigned to the design edges, in matching row order.
#' @param gradient Logical; return the objective gradient instead of its scalar value.
#' @param four Logical; require four rows when TRUE, otherwise at least two.
#' @param positions Strictly increasing variant positions in bases, one per G column.
#' @param x Numeric sample labels to encode as a carrier key.
#' @param alpha Scalar family-wise error probability strictly between zero and one.
#' @param tips Distinct integer tip labels, normally 1:4.
#' @param include_tips Logical; include singleton clades.
#' @param include_root Logical; include the root clade.
#' @param design One ranked-tree design returned by pti_rank_designs.
#' @param rate Positive mutation rate times analyzed length, per generation.
#' @param pair_result Data frame returned by pti_pair_intervals, including its sampling attributes.
#' @param clade Distinct one-based sample indices defining a clade.
#' @returns pti_polymorphic returns variable columns; pti_carrier_keys and pti_key return comma-separated carrier labels. pti_rooted_trees, pti_compatible_trees and pti_tree_clades return lists of nested trees or clades. pti_consensus returns status, compatible_count, supported_clades, trees and a consensus tree. pti_rank_designs returns ranked designs and their edge-duration matrices. pti_fit_rank and pti_fit_ages return likelihood fits, ages, boundaries and optimizer diagnostics. pti_pair_intervals returns pairwise Poisson clock estimates and simultaneous intervals; pti_clade_interval combines these for a clade. pti_linked_constraints returns incompatibility witnesses, minimum_blocks, minimum_changes and cut gaps. pti_poisson_objective returns a scalar negative log likelihood or its gradient; pti_validate_G invisibly returns G.
#' @examples
#' G <- rbind(c(1,1,0), c(1,0,0), c(0,0,1), c(0,0,0))
#' pti_consensus(G)$status
#' pti_pair_intervals(G, mu=1e-4, L=1000)
NULL

#' @rdname models_tree_inference
#' @export
pti_validate_G <- function(G, four = TRUE) {
    if (!is.matrix(G) || !is.numeric(G) || anyNA(G) ||
        any(!is.finite(G)) || any(!G %in% c(0, 1)))
        stop("G must be a finite binary numeric matrix")
    if ((four && nrow(G) != 4L) || (!four && nrow(G) < 2L))
        stop(if (four) "G must contain four sampled copies" else "At least two copies are required")
    invisible(G)
}

#' @rdname models_tree_inference
#' @export
pti_key <- function(x) paste(sort(as.integer(x)), collapse = ",")

#' @rdname models_tree_inference
#' @export
pti_polymorphic <- function(G) {
    pti_validate_G(G)
    keep <- colSums(G) > 0L & colSums(G) < nrow(G)
    G[, keep, drop = FALSE]
}

#' @rdname models_tree_inference
#' @export
pti_rooted_trees <- function(tips = 1:4) {
    if (!is.numeric(tips) || !length(tips) || anyNA(tips) ||
        any(tips != as.integer(tips)) || anyDuplicated(tips))
        stop("tips must be distinct integer labels")
    tips <- sort(as.integer(tips))
    build <- function(x) {
        if (length(x) == 1L) return(list(x))
        result <- list()
        for (mask in 0:(2^(length(x) - 1L) - 2L)) {
            chosen <- as.logical(intToBits(mask)[seq_len(length(x) - 1L)])
            left <- c(x[1], x[-1][chosen])
            right <- setdiff(x, left)
            for (a in build(left)) for (b in build(right))
                result[[length(result) + 1L]] <- list(a, b)
        }
        result
    }
    build(tips)
}

#' @rdname models_tree_inference
#' @export
pti_tree_clades <- function(tree, include_tips = TRUE, include_root = TRUE) {
    recurse <- function(t) {
        if (!is.list(t)) return(list(as.integer(t)))
        c(list(sort(as.integer(unlist(t)))),
          unlist(lapply(t, recurse), recursive = FALSE))
    }
    clades <- recurse(tree)
    root_size <- length(clades[[1]])
    clades[vapply(clades, function(x)
        (include_tips || length(x) > 1L) &&
        (include_root || length(x) < root_size), logical(1))]
}

#' @rdname models_tree_inference
#' @export
pti_carrier_keys <- function(G) {
    if (!ncol(G)) return(character())
    vapply(seq_len(ncol(G)), function(j) pti_key(which(G[, j] == 1)), character(1))
}

#' @rdname models_tree_inference
#' @export
pti_compatible_trees <- function(G) {
    G <- pti_polymorphic(G)
    keys <- unique(pti_carrier_keys(G))
    Filter(function(tree) all(keys %in% vapply(pti_tree_clades(tree), pti_key, character(1))),
           pti_rooted_trees(seq_len(nrow(G))))
}

#' @rdname models_tree_inference
#' @export
pti_consensus <- function(G) {
    G <- pti_polymorphic(G)
    trees <- pti_compatible_trees(G)
    keys <- unique(pti_carrier_keys(G))
    supported <- lapply(keys, function(k) as.integer(strsplit(k, ",", fixed = TRUE)[[1]]))
    supported <- supported[lengths(supported) > 1L & lengths(supported) < nrow(G)]
    supported <- supported[order(lengths(supported), vapply(supported, pti_key, character(1)))]
    if (!length(trees)) return(list(status = "incompatible", compatible_count = 0L,
        supported_clades = supported, trees = list(), tree = NULL))
    clades <- c(lapply(seq_len(nrow(G)), function(i) i), supported, list(seq_len(nrow(G))))
    list(status = if (length(trees) == 1L) "unique_binary" else "unresolved",
         compatible_count = length(trees), supported_clades = supported,
         trees = trees, tree = .cladeTree(clades, seq_len(nrow(G))))
}

.cladeTree <- function(clades, tips) {
    build <- function(parent) {
        if (length(parent) == 1L) return(parent)
        children <- Filter(function(x) length(x) < length(parent) && all(x %in% parent), clades)
        children <- Filter(function(x) !any(vapply(children, function(y)
            length(x) < length(y) && all(x %in% y), logical(1))), children)
        children <- children[order(vapply(children, min, numeric(1)))]
        lapply(children, build)
    }
    build(tips)
}

#' @rdname models_tree_inference
#' @export
pti_rank_designs <- function(tree) {
    tips <- sort(as.integer(unlist(tree)))
    if (!identical(tips, 1:4)) stop("tree must have the four labels 1, 2, 3, 4")
    all_clades <- pti_tree_clades(tree)
    keys <- vapply(all_clades, pti_key, character(1))
    if (length(all_clades) != 7L || anyDuplicated(keys) ||
        any(vapply(Filter(is.list, list(tree)), length, integer(1)) != 2L))
        stop("tree must be rooted and binary")
    # Recursive child counts also exclude a multifurcating nested node.
    binary <- function(t) !is.list(t) || (length(t) == 2L && all(vapply(t, binary, logical(1))))
    if (!binary(tree)) stop("tree must be rooted and binary")
    internal <- all_clades[lengths(all_clades) > 1L]
    internal <- internal[order(lengths(internal), vapply(internal, pti_key, character(1)))]
    ikeys <- vapply(internal, pti_key, character(1))
    root <- pti_key(1:4)
    edge_clades <- all_clades[keys != root]
    edge_clades <- edge_clades[order(lengths(edge_clades), vapply(edge_clades, pti_key, character(1)))]
    ekeys <- vapply(edge_clades, pti_key, character(1))
    parents <- vapply(edge_clades, function(child) {
        possible <- Filter(function(x) length(x) > length(child) && all(child %in% x), all_clades)
        pti_key(possible[[which.min(lengths(possible))]])
    }, character(1))
    perms <- function(x) if (length(x) == 1L) list(x) else
        unlist(lapply(seq_along(x), function(i)
            lapply(perms(x[-i]), function(rest) c(x[i], rest))), recursive = FALSE)
    orders <- Filter(function(order) all(vapply(seq_along(internal), function(i) {
        descendants <- which(vapply(internal, function(x)
            length(x) < length(internal[[i]]) && all(x %in% internal[[i]]), logical(1)))
        all(match(ikeys[descendants], order) < match(ikeys[i], order))
    }, logical(1))), perms(ikeys))
    lapply(orders, function(order) {
        age_row <- function(key) {
            r <- match(key, order)
            if (is.na(r)) return(rep(0, 3L))
            as.numeric(seq_len(3L) <= r)
        }
        A <- t(vapply(seq_along(ekeys), function(i)
            age_row(parents[i]) - age_row(ekeys[i]), numeric(3)))
        dimnames(A) <- list(ekeys, paste0("duration_", 1:3))
        list(tree = tree, order = order, A = A, edge_clades = edge_clades,
             edge_keys = ekeys, parent_keys = parents)
    })
}

#' @rdname models_tree_inference
#' @export
pti_poisson_objective <- function(theta, A, counts, gradient = FALSE) {
    lambda <- as.vector(A %*% theta)
    impossible <- counts > 0 & lambda <= 0
    if (any(impossible) || any(lambda < 0))
        return(if (gradient) rep(-1e100, length(theta)) else 1e100)
    if (gradient) {
        ratio <- numeric(length(counts))
        ratio[counts > 0] <- counts[counts > 0] / lambda[counts > 0]
        return(as.vector(crossprod(A, 1 - ratio)))
    }
    # The 0*log(0) limit is zero; positive counts at zero exposure are impossible.
    positive <- counts > 0
    sum(lambda) - sum(counts[positive] * log(lambda[positive])) + sum(lgamma(counts + 1))
}

#' @rdname models_tree_inference
#' @export
pti_fit_rank <- function(design, counts, rate) {
    A <- design$A
    scale <- max(1, sum(counts) / sum(A))
    solutions <- list()
    # Enumerate faces of the nonnegative orthant, including exact zero edges.
    # The negative log likelihood is convex in theta = mu*L*duration.
    for (mask in 0:7) {
        free <- which(as.logical(intToBits(mask)[1:3]))
        if (any(counts > 0 & rowSums(A[, free, drop = FALSE]) == 0)) next
        if (!length(free)) {
            solutions[[length(solutions) + 1L]] <- list(theta = rep(0, 3),
                objective = pti_poisson_objective(rep(0, 3), A, counts), convergence = 0L)
            next
        }
        starts <- list(rep(scale, length(free)),
                       scale * seq(.35, 1.65, length.out = length(free)),
                       scale * rev(seq(.2, 2, length.out = length(free))))
        expand <- function(x) {theta <- numeric(3); theta[free] <- x; theta}
        for (start in starts) {
            fit <- optim(start, function(x) pti_poisson_objective(expand(x), A, counts),
                gr = function(x) pti_poisson_objective(expand(x), A, counts, TRUE)[free],
                method = "L-BFGS-B", lower = rep(0, length(free)),
                control = list(factr = 1e3, pgtol = 1e-9, maxit = 2000))
            solutions[[length(solutions) + 1L]] <- list(theta = expand(fit$par),
                objective = fit$value, convergence = fit$convergence)
        }
    }
    values <- vapply(solutions, `[[`, numeric(1), "objective")
    best_value <- min(values)
    tied <- which(values <= best_value + 1e-8)
    # Equivalent faces are represented by the solution with exact boundary zeros.
    nz <- vapply(solutions[tied], function(x) sum(x$theta > 0), integer(1))
    best <- solutions[[tied[which.min(nz)]]]
    theta <- best$theta
    gradient <- pti_poisson_objective(theta, A, counts, TRUE)
    positive <- theta > 1e-7
    kkt <- max(c(abs(gradient[positive]), pmax(-gradient[!positive], 0), 0))
    delta <- theta / rate
    ages <- setNames(cumsum(delta), design$order)
    c(design, list(counts = counts, theta = theta, durations = delta, ages = ages,
        branch_lengths = setNames(as.vector(A %*% delta), design$edge_keys),
        expected_counts = setNames(as.vector(A %*% theta), design$edge_keys),
        loglik = -pti_poisson_objective(theta, A, counts), boundary = any(theta == 0),
        kkt_violation = kkt, convergence = best$convergence,
        start_objectives = values))
}

#' @rdname models_tree_inference
#' @export
pti_fit_ages <- function(G, mu, L, tie_tolerance = 1e-7) {
    G <- pti_polymorphic(G)
    if (length(mu) != 1L || !is.finite(mu) || mu <= 0 ||
        length(L) != 1L || !is.finite(L) || L <= 0 ||
        length(tie_tolerance) != 1L || !is.finite(tie_tolerance) || tie_tolerance < 0)
        stop("mu, L, and tie_tolerance must be finite; mu and L must be positive")
    consensus <- pti_consensus(G)
    if (!length(consensus$trees)) return(list(status = "incompatible", consensus = consensus,
        fits = list(), best_fits = list(), node_ages = data.frame()))
    observed_keys <- pti_carrier_keys(G)
    fits <- list()
    for (tree_index in seq_along(consensus$trees)) {
        for (design in pti_rank_designs(consensus$trees[[tree_index]])) {
            counts <- vapply(design$edge_keys, function(k) sum(observed_keys == k), numeric(1))
            fit <- pti_fit_rank(design, counts, mu * L)
            fit$tree_index <- tree_index
            fits[[length(fits) + 1L]] <- fit
        }
    }
    best_loglik <- max(vapply(fits, `[[`, numeric(1), "loglik"))
    best_fits <- Filter(function(x) x$loglik >= best_loglik - tie_tolerance, fits)
    reported <- c(vapply(consensus$supported_clades, pti_key, character(1)), pti_key(1:4))
    node_ages <- do.call(rbind, lapply(reported, function(key) {
        ages <- vapply(best_fits, function(x) unname(x$ages[key]), numeric(1))
        data.frame(clade = key, estimate_min = min(ages), estimate_max = max(ages),
            same_height_across_best_fits = diff(range(ages)) <= 1e-5 * max(1, max(ages)),
            row.names = NULL)
    }))
    list(status = consensus$status, consensus = consensus, fits = fits,
         best_fits = best_fits, best_loglik = best_loglik, node_ages = node_ages,
         mu = mu, analyzed_length_bp = L, assumption = "one error-free infinite-sites genealogy")
}

#' @rdname models_tree_inference
#' @export
pti_pair_intervals <- function(G, mu, L, alpha = .05) {
    pti_validate_G(G, four = FALSE)
    if (length(mu) != 1L || !is.finite(mu) || mu <= 0 ||
        length(L) != 1L || !is.finite(L) || L <= 0 ||
        length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
        stop("positive finite mu and L and 0 < alpha < 1 are required")
    pairs <- combn(seq_len(nrow(G)), 2)
    result <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
        i <- pairs[1, k]; j <- pairs[2, k]
        count <- sum(G[i, ] != G[j, ])
        ci <- poisson.test(count, T = 2 * mu * L,
                           conf.level = 1 - alpha / ncol(pairs))$conf.int
        data.frame(i = i, j = j, differences = count, estimate = count / (2 * mu * L),
                   lower = ci[1], upper = ci[2])
    }))
    attr(result, "alpha") <- alpha
    attr(result, "pair_count") <- ncol(pairs)
    attr(result, "sample_count") <- nrow(G)
    attr(result, "assumption") <- "one error-free infinite-sites genealogy; known mutation rate"
    result
}

#' @rdname models_tree_inference
#' @export
pti_clade_interval <- function(pair_result, clade = 1:4) {
    if (!all(c("i", "j", "estimate", "lower", "upper") %in% names(pair_result)))
        stop("pair_result must come from pti_pair_intervals")
    if (!length(clade) || anyNA(clade) || anyDuplicated(clade) ||
        any(!clade %in% seq_len(attr(pair_result, "sample_count")))) stop("invalid clade")
    if (length(clade) == 1L) return(list(clade = pti_key(clade), estimate_max_pair = 0,
        lower = 0, upper = 0, simultaneous_coverage_at_least = 1 - attr(pair_result, "alpha")))
    x <- pair_result[pair_result$i %in% clade & pair_result$j %in% clade, , drop = FALSE]
    list(clade = pti_key(clade), estimate_max_pair = max(x$estimate), lower = max(x$lower),
         upper = max(x$upper), simultaneous_coverage_at_least = 1 - attr(pair_result, "alpha"))
}

#' @rdname models_tree_inference
#' @export
pti_linked_constraints <- function(G, positions) {
    pti_validate_G(G)
    if (!is.numeric(positions) || length(positions) != ncol(G) || anyNA(positions) ||
        any(!is.finite(positions)) || is.unsorted(positions, strictly = TRUE))
        stop("positions must be finite, strictly increasing, and match the columns of G")
    keep <- which(colSums(G) > 0L & colSums(G) < nrow(G))
    G <- G[, keep, drop = FALSE]; positions <- positions[keep]
    m <- ncol(G)
    keys <- pti_carrier_keys(G)
    sets <- lapply(seq_len(m), function(j) which(G[, j] == 1))
    incompatible <- matrix(FALSE, m, m)
    witnesses <- list()
    if (m >= 2L) for (a in seq_len(m - 1L)) for (b in seq.int(a + 1L, m)) {
        A <- sets[[a]]; B <- sets[[b]]
        conflict <- length(intersect(A, B)) > 0L && length(setdiff(A, B)) > 0L &&
                    length(setdiff(B, A)) > 0L
        incompatible[a, b] <- incompatible[b, a] <- conflict
        if (conflict) witnesses[[length(witnesses) + 1L]] <- data.frame(
            left_site = keep[a], right_site = keep[b], left = positions[a], right = positions[b],
            left_clade = keys[a], right_clade = keys[b])
    }
    witnesses <- if (length(witnesses)) do.call(rbind, witnesses) else
        data.frame(left_site = integer(), right_site = integer(), left = numeric(),
                   right = numeric(), left_clade = character(), right_clade = character())
    disjoint <- witnesses[FALSE, , drop = FALSE]
    if (nrow(witnesses)) {
        candidates <- witnesses[order(witnesses$right, -witnesses$left), , drop = FALSE]
        last_right <- -Inf
        for (i in seq_len(nrow(candidates))) if (candidates$left[i] >= last_right) {
            disjoint <- rbind(disjoint, candidates[i, , drop = FALSE])
            last_right <- candidates$right[i]
        }
    }
    if (!m) return(list(witnesses = witnesses, disjoint_witnesses = disjoint,
        minimum_blocks = 1L, minimum_changes = 0L,
        cut_gaps = data.frame(after_site = integer(), before_site = integer(), left = numeric(), right = numeric()),
        minimal_cut_sets = list(integer()), informative_sites = keep,
        compatible_blocks = matrix(logical(), 0, 0), positions = positions))
    compatible <- matrix(FALSE, m, m)
    for (i in seq_len(m)) for (j in seq.int(i, m))
        compatible[i, j] <- !any(incompatible[i:j, i:j, drop = FALSE])
    prefix <- c(0, rep(Inf, m))
    previous <- vector("list", m)
    for (j in seq_len(m)) {
        starts <- which(compatible[, j])
        scores <- prefix[starts] + 1
        prefix[j + 1L] <- min(scores)
        previous[[j]] <- starts[scores == min(scores)] - 1L
    }
    recover <- function(j) {
        if (!j) return(list(integer()))
        unlist(lapply(previous[[j]], function(cut)
            lapply(recover(cut), function(rest) if (cut) c(rest, cut) else rest)), recursive = FALSE)
    }
    minimal_cut_sets <- recover(m)
    all_cuts <- sort(unique(unlist(minimal_cut_sets)))
    cut_gaps <- if (length(all_cuts)) data.frame(after_site = keep[all_cuts],
        before_site = keep[all_cuts + 1L], left = positions[all_cuts], right = positions[all_cuts + 1L]) else
        data.frame(after_site = integer(), before_site = integer(), left = numeric(), right = numeric())
    list(witnesses = witnesses, disjoint_witnesses = disjoint,
         minimum_blocks = as.integer(prefix[m + 1L]),
         minimum_changes = as.integer(prefix[m + 1L] - 1L), cut_gaps = cut_gaps,
         minimal_cut_sets = minimal_cut_sets, informative_sites = keep,
         compatible_blocks = compatible, positions = positions,
         interpretation = "minimum piecewise-compatible explanation; cut gaps are not confidence intervals")
}

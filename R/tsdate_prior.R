# From Daniel Nachun's tskitR: .compute_conditional_prior,
# .spans_by_samples and build_prior_grid (tsdate-prior.R).
# GPL >= 3.

tdLognormApprox <- function(mean, var) {
    beta <- log(var/mean^2 + 1)
    alpha <- log(mean) - 0.5 * beta
    c(alpha, beta)
}

tdGammaApprox <- function(mean, variance) {
    c(mean^2/variance, mean/variance)
}

tdTauExpect <- function(i, n) {
    ifelse(i == n, 2 * (1 - 1/n), (i - 1)/n)
}

tdTauVarMrca <- function(n) {
    v <- seq(2, n)
    var_val <- sum(1/(v^2 * (v - 1)^2))
    abs(4 * var_val)
}

tdComputeConditionalPrior <- function(total_tips, prior_distr = "lognorm", approx_priors = NULL) {
    func_approx <- if (prior_distr == "lognorm") 
        tdLognormApprox
    else tdGammaApprox
    priors <- matrix(NA_real_, nrow = total_tips + 1L, ncol = 4L)
    colnames(priors) <- c("alpha", "beta", "mean", "var")
    if (prior_distr == "lognorm") {
        priors[2, ] <- c(-Inf, 0, 0, 0)
    }
    else {
        priors[2, ] <- c(0, 1, 0, 0)
    }
    all_tips <- seq(2L, total_tips)
    if (!is.null(approx_priors)) {
        variances <- approx(approx_priors[, 1], approx_priors[, 2], xout = all_tips/total_tips, rule = 2)$y
        variances[all_tips == total_tips] <- tdTauVarMrca(total_tips)
    }
    else {
        var_all <- tdTsdateConditionalCoalescentVariance(total_tips)
        variances <- var_all[all_tips + 1L]
    }
    for (idx in seq_along(all_tips)) {
        tips <- all_tips[idx]
        expectation <- tdTauExpect(tips, total_tips)
        var_val <- variances[idx]
        if (is.na(var_val) || var_val <= 0 || expectation <= 0) {
            next
        }
        params <- func_approx(expectation, var_val)
        priors[tips + 1L, ] <- c(params[1], params[2], expectation, var_val)
    }
    priors
}

tdSpansBySamples <- function(ts) {
    samples <- tdSamples(ts)
    num_nodes <- tdNumNodes(ts)
    sample_set <- samples
    spans <- list()
    node_spans <- numeric(num_nodes)
    tree_obj <- tdTree(ts)
    tdFirst(tree_obj)
    n_tips_per_tree <- integer(tdNumTrees(ts))
    tree_idx <- 1L
    repeat {
        interval <- tdInterval(tree_obj)
        span <- interval[2] - interval[1]
        n_samples <- 0L
        roots <- tdRoots(tree_obj)
        for (s in sample_set) {
            if (tdParentNode(tree_obj, s) != -1L || s %in% roots) {
                n_samples <- n_samples + 1L
            }
        }
        n_tips_per_tree[tree_idx] <- n_samples
        all_nodes <- tdPreorder(tree_obj)
        for (node in all_nodes) {
            if (node %in% sample_set) {
                node_spans[node + 1L] <- node_spans[node + 1L] + span
                next
            }
            n_desc <- tdNumSamplesNode(tree_obj, node)
            if (n_desc > 0L) {
                node_spans[node + 1L] <- node_spans[node + 1L] + span
                node_key <- as.character(node)
                tips_key <- as.character(n_samples)
                if (is.null(spans[[node_key]])) {
                  spans[[node_key]] <- list()
                }
                if (is.null(spans[[node_key]][[tips_key]])) {
                  spans[[node_key]][[tips_key]] <- list()
                }
                desc_key <- as.character(n_desc)
                prev <- spans[[node_key]][[tips_key]][[desc_key]]
                if (is.null(prev)) 
                  prev <- 0
                spans[[node_key]][[tips_key]][[desc_key]] <- prev + span
            }
        }
        tree_idx <- tree_idx + 1L
        if (!tdNextTree(tree_obj)) 
            break
    }
    nodes_to_date <- as.integer(names(spans))
    total_fixed_at_0_counts <- unique(n_tips_per_tree)
    total_fixed_at_0_counts <- total_fixed_at_0_counts[total_fixed_at_0_counts > 0]
    node_span_data <- list()
    for (node_key in names(spans)) {
        node_span_data[[node_key]] <- list()
        for (tips_key in names(spans[[node_key]])) {
            entries <- spans[[node_key]][[tips_key]]
            desc_tips <- as.integer(names(entries))
            span_vals <- as.numeric(unlist(entries))
            node_span_data[[node_key]][[tips_key]] <- data.frame(descendant_tips = desc_tips, span = span_vals, 
                stringsAsFactors = FALSE)
        }
    }
    list(node_span_data = node_span_data, nodes_to_date = nodes_to_date, node_spans = node_spans, total_fixed_at_0_counts = total_fixed_at_0_counts)
}

tdMixtureExpectAndVar <- function(mixture, prior_store) {
    expectation <- 0
    first <- 0
    secnd <- 0
    weight_sum <- 0
    for (total_tips_key in names(mixture)) {
        tip_df <- mixture[[total_tips_key]]
        total_tips <- as.integer(total_tips_key)
        prior_mat <- prior_store[[total_tips_key]]
        desc_indices <- tip_df$descendant_tips + 1L
        mean_vals <- prior_mat[desc_indices, 3]
        var_vals <- prior_mat[desc_indices, 4]
        w <- tip_df$span
        expectation <- expectation + sum(mean_vals * w, na.rm = TRUE)
        first <- first + sum(var_vals * w, na.rm = TRUE)
        secnd <- secnd + sum(mean_vals^2 * w, na.rm = TRUE)
        weight_sum <- weight_sum + sum(w)
    }
    mean_out <- expectation/weight_sum
    var_out <- (first + secnd)/weight_sum - mean_out^2
    c(mean_out, var_out)
}

tdCreateTimepoints <- function(prior_params, n_points = 21L, prior_distr = "lognorm") {
    if (prior_distr == "lognorm") {
        ppf <- function(p, alpha, beta) {
            qlnorm(p, meanlog = alpha, sdlog = sqrt(beta))
        }
        cdf <- function(t, alpha, beta) {
            plnorm(t, meanlog = alpha, sdlog = sqrt(beta))
        }
    }
    else {
        ppf <- function(p, alpha, beta) {
            qgamma(p, shape = alpha, rate = beta)
        }
        cdf <- function(t, alpha, beta) {
            pgamma(t, shape = alpha, rate = beta)
        }
    }
    percentiles <- seq(0, 1, length.out = n_points + 1L)[2:n_points]
    alpha0 <- prior_params[3, 1]
    beta0 <- prior_params[3, 2]
    t_set <- ppf(percentiles, alpha0, beta0)
    max_tips <- nrow(prior_params)
    max_sep <- 1/(n_points - 1L)
    if (max_tips > 3) {
        for (i in seq(4, max_tips)) {
            alpha_i <- prior_params[i, 1]
            beta_i <- prior_params[i, 2]
            if (is.na(alpha_i) || is.na(beta_i)) 
                next
            proj <- cdf(t_set, alpha_i, beta_i)
            dists <- vapply(percentiles, function(v) min(abs(v - proj)), numeric(1))
            wd <- which(dists > max_sep)
            if (length(wd) > 0) {
                new_points <- ppf(percentiles[wd], alpha_i, beta_i)
                t_set <- c(t_set, new_points)
            }
        }
    }
    t_set <- sort(t_set)
    c(0, t_set)
}

tdFillPriors <- function(prior_params, timepoints, ts, pop_hist, prior_distr = "lognorm") {
    samps <- tdSamples(ts)
    n_nodes <- tdNumNodes(ts)
    datable_nodes <- setdiff(seq(0L, n_nodes - 1L), samps)
    tc <- tdTables(ts)
    nt <- tdNodes(tc)
    node_times <- nt$time
    datable_sorted <- datable_nodes[order(node_times[datable_nodes + 1L])]
    tp_natural <- tdToNaturalTimescale(pop_hist, timepoints)
    ntv <- tdNodeTimeValues(n_nodes, as.integer(datable_sorted), tp_natural)
    for (node in datable_sorted) {
        alpha <- prior_params[node + 1L, 1]
        beta <- prior_params[node + 1L, 2]
        if (is.na(alpha) || is.na(beta)) 
            next
        if (prior_distr == "lognorm") {
            prior_node <- plnorm(timepoints, meanlog = alpha, sdlog = sqrt(beta))
        }
        else {
            prior_node <- pgamma(timepoints, shape = alpha, scale = 1/beta)
        }
        max_val <- max(prior_node, na.rm = TRUE)
        if (max_val > 0) {
            prior_node <- prior_node/max_val
        }
        prior_diff <- c(0, diff(prior_node))
        ntv <- tdNtvSet(ntv, node, prior_diff)
    }
    ntv <- tdNtvStandardize(ntv)
    ntv
}

tdBuildPriorGrid <- function(ts, population_size, timepoints = 20L, approximate_priors = FALSE, approx_prior_size = NULL, 
    prior_distribution = "lognorm") {
    if (is.numeric(population_size) && !tdIs(population_size, "PopulationSizeHistory")) {
        pop_hist <- tdPopulationSizeHistory(population_size)
    }
    else {
        pop_hist <- population_size
    }
    span_data <- tdSpansBySamples(ts)
    total_counts <- span_data$total_fixed_at_0_counts
    n_samples <- length(tdSamples(ts))
    all_counts <- unique(c(n_samples, total_counts))
    all_counts <- all_counts[all_counts > 0]
    prior_store <- list()
    for (n in all_counts) {
        approx_table <- if (approximate_priors && !is.null(approx_prior_size)) {
            var_exact <- tdTsdateConditionalCoalescentVariance(approx_prior_size + 1L)
            all_t <- seq(2, approx_prior_size)
            cbind(all_t/approx_prior_size, var_exact[all_t + 1L])
        }
        else {
            NULL
        }
        prior_store[[as.character(n)]] <- tdComputeConditionalPrior(n, prior_distribution, approx_table)
    }
    func_approx <- if (prior_distribution == "lognorm") 
        tdLognormApprox
    else tdGammaApprox
    n_nodes <- tdNumNodes(ts)
    prior_params <- matrix(NA_real_, nrow = n_nodes, ncol = 4L)
    colnames(prior_params) <- c("alpha", "beta", "mean", "var")
    for (node_key in names(span_data$node_span_data)) {
        node <- as.integer(node_key)
        mixture <- span_data$node_span_data[[node_key]]
        mv <- tdMixtureExpectAndVar(mixture, prior_store)
        if (mv[1] > 0 && mv[2] > 0) {
            params <- func_approx(mv[1], mv[2])
            prior_params[node + 1L, ] <- c(params[1], params[2], mv[1], mv[2])
        }
    }
    if (is.numeric(timepoints) && length(timepoints) == 1 && timepoints == as.integer(timepoints)) {
        max_n <- max(all_counts)
        base_prior <- prior_store[[as.character(max_n)]]
        tp <- tdCreateTimepoints(base_prior, as.integer(timepoints) + 1L, prior_distribution)
    }
    else {
        tp <- sort(as.numeric(timepoints))
        tp <- tdToCoalescentTimescale(pop_hist, tp)
    }
    tdFillPriors(prior_params, tp, ts, pop_hist, prior_distribution)
}

# Conditional coalescent moments from the supplied tsdate_prior.cpp.

tdMarginalizeOverAncestors <- function(val) {
    n <- nrow(val)
    out <- matrix(0, n + 1L, ncol(val))
    logProbability <- c(NaN, NaN, 0)
    if (n >= 3L) for (k in seq.int(n - 1L, 2L)) {
        constant <- log(n - k) + log(k - 2) - log(k + 1)
        for (a in 2:(n - k + 1L)) {
            out[k + 1L, ] <- out[k + 1L, ] + exp(logProbability[a + 1L]) * val[a + 1L, ]
            if (k > 2L) logProbability[a + 1L] <-
                logProbability[a + 1L] + constant - log(n - a - k + 2)
        }
        if (k > 2L) logProbability <- c(logProbability,
            tail(logProbability, 1L) + log(n - k + 2) - log(k + 1) - constant)
    }
    out[n + 1L, ] <- val[2L, ]
    out
}
tdTsdateConditionalCoalescentVariance <- function(num_tips) {
    rates <- c(0, 2 / ((2:num_tips) * (1:(num_tips - 1L))))
    # Hypoexponential moments sum every wait, including the final array entry.
    # The supplied C++ reverse loop skips the penultimate entry. This boundary
    # correction is tested against the analytical MRCA variance, including 10/9
    # for three tips, and the source tsdate cumulative-sum definition.
    mean <- rev(cumsum(rev(rates)))
    variance <- rev(cumsum(rev(rates^2)))
    moments <- tdMarginalizeOverAncestors(cbind(mean, variance + mean^2))
    moments[, 2L] - moments[, 1L]^2
}

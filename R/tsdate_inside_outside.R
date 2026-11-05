# Source functions copied from Daniel Nachun's tskitR/R/tsdate-discrete.R.
# GPL >= 3. Exact source hashes and function names are in source-map.tsv.
# Mechanical adaptation: private .td_ names and R-list table accessors.
# The inside/outside recurrences and conditional coalescent priors are retained.

.td_LIN_GRID <- "linear"

.td_LOG_GRID <- "logarithmic"

.td_GAMMA_PAR <- "gamma_parameter"

.td_node_time_values <- function(num_nodes, nonfixed_nodes, timepoints, fill_value = NaN) {
    stopifnot(is.integer(nonfixed_nodes) || is.numeric(nonfixed_nodes))
    nonfixed_nodes <- as.integer(nonfixed_nodes)
    stopifnot(all(nonfixed_nodes >= 0L & nonfixed_nodes < num_nodes))
    if (is.numeric(timepoints) && length(timepoints) > 1) {
        grid_size <- length(timepoints)
        probability_space <- .td_LIN_GRID
    }
    else if (length(timepoints) == 1 && timepoints == as.integer(timepoints)) {
        grid_size <- as.integer(timepoints)
        probability_space <- .td_GAMMA_PAR
    }
    else {
        grid_size <- length(timepoints)
        probability_space <- .td_LIN_GRID
    }
    num_nonfixed <- length(nonfixed_nodes)
    num_fixed <- num_nodes - num_nonfixed
    grid_data <- matrix(fill_value, nrow = num_nonfixed, ncol = grid_size)
    fixed_data <- rep(fill_value, num_fixed)
    row_lookup <- integer(num_nodes)
    row_lookup[nonfixed_nodes + 1L] <- seq_along(nonfixed_nodes) - 1L
    fixed_nodes <- setdiff(seq(0L, num_nodes - 1L), nonfixed_nodes)
    if (length(fixed_nodes) > 0) {
        row_lookup[fixed_nodes + 1L] <- -(seq_along(fixed_nodes))
    }
    obj <- list(num_nodes = num_nodes, nonfixed_nodes = nonfixed_nodes, num_nonfixed = num_nonfixed, 
        grid_data = grid_data, fixed_data = fixed_data, row_lookup = row_lookup, timepoints = timepoints, 
        probability_space = probability_space, grid_size = grid_size)
    class(obj) <- "NodeTimeValues"
    obj
}

.td_ntv_get <- function(ntv, node_id) {
    idx <- ntv$row_lookup[node_id + 1L]
    if (idx < 0) {
        ntv$fixed_data[-idx]
    }
    else {
        ntv$grid_data[idx + 1L, ]
    }
}

.td_ntv_set <- function(ntv, node_id, value) {
    idx <- ntv$row_lookup[node_id + 1L]
    if (idx < 0) {
        ntv$fixed_data[-idx] <- value
    }
    else {
        ntv$grid_data[idx + 1L, ] <- value
    }
    ntv
}

.td_ntv_is_fixed <- function(ntv, node_id) {
    ntv$row_lookup[node_id + 1L] < 0
}

.td_ntv_clone <- function(ntv, grid_data = NaN, fixed_data = NULL, probability_space = NULL) {
    new_ntv <- list(num_nodes = ntv$num_nodes, nonfixed_nodes = ntv$nonfixed_nodes, num_nonfixed = ntv$num_nonfixed, 
        row_lookup = ntv$row_lookup, timepoints = ntv$timepoints, grid_size = ntv$grid_size)
    if (is.matrix(grid_data)) {
        stopifnot(identical(dim(grid_data), dim(ntv$grid_data)))
        new_ntv$grid_data <- grid_data
        if (is.null(fixed_data)) {
            new_ntv$fixed_data <- rep(NaN, length(ntv$fixed_data))
        }
        else if (length(fixed_data) == 1) {
            new_ntv$fixed_data <- rep(fixed_data, length(ntv$fixed_data))
        }
        else {
            new_ntv$fixed_data <- fixed_data
        }
    }
    else {
        new_ntv$grid_data <- matrix(grid_data, nrow = nrow(ntv$grid_data), ncol = ncol(ntv$grid_data))
        if (is.null(fixed_data)) {
            new_ntv$fixed_data <- rep(grid_data, length(ntv$fixed_data))
        }
        else if (length(fixed_data) == 1) {
            new_ntv$fixed_data <- rep(fixed_data, length(ntv$fixed_data))
        }
        else {
            new_ntv$fixed_data <- fixed_data
        }
    }
    if (is.null(probability_space)) {
        new_ntv$probability_space <- ntv$probability_space
    }
    else {
        new_ntv$probability_space <- probability_space
    }
    class(new_ntv) <- "NodeTimeValues"
    new_ntv
}

.td_ntv_force_probability_space <- function(ntv, target_space) {
    if (target_space == .td_LIN_GRID) {
        if (ntv$probability_space == .td_LIN_GRID) {
            return(ntv)
        }
        else if (ntv$probability_space == .td_LOG_GRID) {
            ntv$grid_data <- exp(ntv$grid_data)
            ntv$fixed_data <- exp(ntv$fixed_data)
            ntv$probability_space <- .td_LIN_GRID
            return(ntv)
        }
    }
    else if (target_space == .td_LOG_GRID) {
        if (ntv$probability_space == .td_LOG_GRID) {
            return(ntv)
        }
        else if (ntv$probability_space == .td_LIN_GRID) {
            ntv$grid_data <- log(ntv$grid_data)
            ntv$fixed_data <- log(ntv$fixed_data)
            ntv$probability_space <- .td_LOG_GRID
            return(ntv)
        }
    }
    else if (target_space == .td_GAMMA_PAR) {
        if (ntv$probability_space == .td_GAMMA_PAR) {
            return(ntv)
        }
    }
    stop(sprintf("Cannot force %s probabilities into %s space", ntv$probability_space, target_space))
}

.td_ntv_standardize <- function(ntv) {
    if (ncol(ntv$grid_data) > 1) {
        rowmax <- apply(ntv$grid_data[, -1, drop = FALSE], 1, max)
    }
    else {
        rowmax <- ntv$grid_data[, 1]
    }
    if (ntv$probability_space == .td_LIN_GRID) {
        ntv$grid_data <- ntv$grid_data/rowmax
    }
    else if (ntv$probability_space == .td_LOG_GRID) {
        ntv$grid_data <- ntv$grid_data - rowmax
    }
    else {
        stop("Cannot standardize in probability space: ", ntv$probability_space)
    }
    ntv
}

.td_ntv_to_probabilities <- function(ntv) {
    stopifnot(ntv$probability_space == .td_LIN_GRID)
    row_sums <- rowSums(ntv$grid_data)
    ntv$grid_data <- ntv$grid_data/row_sums
    ntv
}

.td_discrete_likelihoods <- function(ts, timepoints, mutation_rate = NULL, eps = 0, fixed_node_set = NULL, 
    standardize = FALSE, log_space = FALSE) {
    stopifnot(.td_is(ts, "TreeSequence"))
    timepoints <- as.numeric(timepoints)
    grid_size <- length(timepoints)
    if (is.null(fixed_node_set)) {
        fixed_node_set <- as.integer(.td_samples(ts))
    }
    fixed_set <- as.integer(fixed_node_set)
    mut_edges <- .td_count_edge_mutations(ts)
    if (log_space) {
        identity_constant <- 0
        null_constant <- -Inf
        probability_space <- .td_LOG_GRID
    }
    else {
        identity_constant <- 1
        null_constant <- 0
        probability_space <- .td_LIN_GRID
    }
    timediff_lower_tri <- unlist(lapply(seq_along(timepoints), function(i) {
        timepoints[i] - timepoints[1:i] + eps
    }))
    timediff <- timepoints - timepoints[1] + eps
    tri_size <- as.integer(grid_size * (grid_size + 1)/2)
    row_indices_starts <- integer(grid_size)
    for (t in seq_len(grid_size)) {
        row_indices_starts[t] <- as.integer((t - 1L) * t/2)
    }
    col_indices <- integer(grid_size)
    running_sum <- 0L
    for (i in seq_len(grid_size)) {
        col_indices[i] <- running_sum
        running_sum <- running_sum + grid_size - i + 1L
    }
    to_lower_tri <- unlist(lapply(seq_len(grid_size), function(i) seq_len(i)))
    to_upper_tri <- unlist(lapply(seq_len(grid_size), function(i) seq(i, grid_size)))
    row_indices_for_reorder <- lapply(seq_len(grid_size), function(t) {
        nn <- seq(t, grid_size)
        as.integer((nn * (nn - 1))/2 + t)
    })
    reorder_indices <- unlist(row_indices_for_reorder)
    tc <- .td_tables(ts)
    et <- .td_edges(tc)
    edge_spans <- et$right - et$left
    edge_children <- et$child
    num_edges <- .td_num_rows(et)
    obj <- list(ts = ts, timepoints = timepoints, grid_size = grid_size, tri_size = tri_size, fixed_set = fixed_set, 
        mut_rate = mutation_rate, mut_edges = mut_edges, standardize = standardize, log_space = log_space, 
        identity_constant = identity_constant, null_constant = null_constant, probability_space = probability_space, 
        timediff_lower_tri = timediff_lower_tri, timediff = timediff, row_indices_starts = row_indices_starts, 
        col_indices = col_indices, to_lower_tri = to_lower_tri, to_upper_tri = to_upper_tri, reorder_indices = reorder_indices, 
        edge_spans = edge_spans, edge_children = edge_children, num_edges = num_edges, unfixed_lik_cache = NULL)
    if (log_space) {
        class(obj) <- c("LogLikelihoods", "Likelihoods")
    }
    else {
        class(obj) <- "Likelihoods"
    }
    obj
}

.td_lik_precalculate <- function(lik) {
    if (is.null(lik$mut_rate)) {
        stop("Cannot calculate mutation likelihoods with no mutation_rate set")
    }
    is_fixed_child <- lik$edge_children %in% lik$fixed_set
    result <- .td_tsdate_precalc_mutation_liks(as.integer(lik$mut_edges), as.numeric(lik$edge_spans), 
        as.logical(is_fixed_child), as.numeric(lik$timediff_lower_tri), as.numeric(lik$mut_rate), as.logical(lik$log_space), 
        as.logical(lik$standardize))
    lik$unfixed_lik_cache <- result$cache
    lik$edge_cache_idx <- result$edge_cache_idx
    lik
}

.td_lik_get_mut_lower_tri <- function(lik, edge_id) {
    cache_idx <- lik$edge_cache_idx[edge_id + 1L]
    stopifnot(cache_idx >= 0)
    lik$unfixed_lik_cache[cache_idx + 1L, ]
}

.td_lik_get_mut_upper_tri <- function(lik, edge_id) {
    lower <- .td_lik_get_mut_lower_tri(lik, edge_id)
    lower[lik$reorder_indices]
}

.td_lik_get_mut_fixed <- function(lik, edge_id) {
    muts <- lik$mut_edges[edge_id + 1L]
    span <- lik$edge_spans[edge_id + 1L]
    .td_tsdate_mutation_lik_fixed(as.integer(muts), as.numeric(span), as.numeric(lik$timediff), as.numeric(lik$mut_rate), 
        as.logical(lik$log_space), as.logical(lik$standardize))
}

.td_lik_make_lower_tri <- function(lik, arr) {
    arr[lik$to_lower_tri]
}

.td_lik_make_upper_tri <- function(lik, arr) {
    arr[lik$to_upper_tri]
}

.td_lik_rowsum_lower_tri <- function(lik, arr) {
    if (lik$log_space) {
        .td_tsdate_rowsum_lower_tri_log(as.numeric(arr), as.integer(lik$row_indices_starts))
    }
    else {
        .td_tsdate_rowsum_lower_tri_lin(as.numeric(arr), as.integer(lik$row_indices_starts))
    }
}

.td_lik_rowsum_upper_tri <- function(lik, arr) {
    if (lik$log_space) {
        .td_tsdate_rowsum_upper_tri_log(as.numeric(arr), as.integer(lik$col_indices))
    }
    else {
        .td_tsdate_rowsum_upper_tri_lin(as.numeric(arr), as.integer(lik$col_indices))
    }
}

.td_lik_combine <- function(lik, a, b) {
    if (lik$log_space) 
        a + b
    else a * b
}

.td_lik_ratio <- function(lik, a, b, div_0_null = FALSE) {
    if (lik$log_space) {
        ret <- a - b
    }
    else {
        ret <- a/b
    }
    if (div_0_null) {
        ret[is.nan(ret)] <- lik$null_constant
    }
    ret
}

.td_lik_marginalize <- function(lik, arr) {
    if (lik$log_space) {
        .td_tsdate_logsumexp(as.numeric(arr))
    }
    else {
        sum(arr)
    }
}

.td_lik_scale_geometric <- function(lik, fraction, value) {
    if (lik$log_space) {
        fraction * value
    }
    else {
        value^fraction
    }
}

.td_lik_get_inside <- function(lik, arr, edge_id) {
    liks <- lik$identity_constant
    if (!is.null(lik$mut_rate)) {
        mut_lik <- .td_lik_get_mut_lower_tri(lik, edge_id)
        liks <- .td_lik_combine(lik, liks, mut_lik)
    }
    .td_lik_rowsum_lower_tri(lik, .td_lik_combine(lik, arr, liks))
}

.td_lik_get_outside <- function(lik, arr, edge_id) {
    liks <- lik$identity_constant
    if (!is.null(lik$mut_rate)) {
        mut_lik <- .td_lik_get_mut_upper_tri(lik, edge_id)
        liks <- .td_lik_combine(lik, liks, mut_lik)
    }
    .td_lik_rowsum_upper_tri(lik, .td_lik_combine(lik, arr, liks))
}

.td_lik_get_fixed <- function(lik, arr, edge_id) {
    liks <- lik$identity_constant
    if (!is.null(lik$mut_rate)) {
        mut_lik <- .td_lik_get_mut_fixed(lik, edge_id)
        liks <- .td_lik_combine(lik, liks, mut_lik)
    }
    .td_lik_combine(lik, arr, liks)
}

.td_belief_propagation <- function(priors, lik) {
    nonfixed <- priors$nonfixed_nodes
    fixed <- lik$fixed_set
    overlap <- intersect(nonfixed, fixed)
    if (length(overlap) > 0) {
        stop("Prior and likelihood objects disagree on which nodes are fixed")
    }
    if (length(fixed) + length(nonfixed) != .td_num_nodes(lik$ts)) {
        stop("Prior and likelihood objects disagree on which nodes are fixed")
    }
    priors <- .td_ntv_force_probability_space(priors, lik$probability_space)
    tc <- .td_tables(lik$ts)
    et <- .td_edges(tc)
    edge_children <- et$child
    edge_parents <- et$parent
    edge_lefts <- et$left
    edge_rights <- et$right
    edge_ids <- seq(0L, .td_num_rows(et) - 1L)
    edge_spans <- edge_rights - edge_lefts
    num_nodes <- .td_num_nodes(lik$ts)
    spans <- numeric(num_nodes)
    for (i in seq_along(edge_children)) {
        spans[edge_children[i] + 1L] <- spans[edge_children[i] + 1L] + edge_spans[i]
    }
    root_spans <- list()
    tr <- .td_tree(lik$ts)
    .td_first(tr)
    repeat {
        iv <- .td_interval(tr)
        span <- iv[2] - iv[1]
        r <- .td_roots(tr)
        if (length(r) == 1L) {
            key <- as.character(r)
            prev <- root_spans[[key]]
            if (is.null(prev)) 
                prev <- 0
            root_spans[[key]] <- prev + span
        }
        if (!.td_next_tree(tr)) 
            break
    }
    for (key in names(root_spans)) {
        node <- as.integer(key)
        spans[node + 1L] <- spans[node + 1L] + root_spans[[key]]
    }
    obj <- list(priors = priors, lik = lik, ts = lik$ts, fixed_set = lik$fixed_set, nonfixed_nodes = nonfixed, 
        spans = spans, root_spans = root_spans, edge_children = edge_children, edge_parents = edge_parents, 
        edge_lefts = edge_lefts, edge_rights = edge_rights, edge_ids = edge_ids, edge_spans = edge_spans, 
        num_edges = length(edge_children), inside = NULL, outside = NULL, denominator = NULL, posterior_grid = NULL, 
        posterior_mean = NULL)
    class(obj) <- "BeliefPropagation"
    obj
}

.td_edges_by_parent_asc <- function(bp) {
    parents <- bp$edge_parents
    unique_parents <- unique(parents)
    groups <- list()
    for (p in unique_parents) {
        mask <- which(parents == p)
        groups[[as.character(p)]] <- list(parent = p, edge_indices = mask - 1L)
    }
    tc <- .td_tables(bp$ts)
    nt <- .td_nodes(tc)
    node_times <- nt$time
    parent_times <- node_times[unique_parents + 1L]
    ord <- order(parent_times)
    groups[as.character(unique_parents[ord])]
}

.td_edges_by_child_desc <- function(bp) {
    children <- bp$edge_children
    unique_children <- unique(children)
    tc <- .td_tables(bp$ts)
    nt <- .td_nodes(tc)
    node_times <- nt$time
    groups <- list()
    for (c_node in unique_children) {
        mask <- which(children == c_node)
        groups[[as.character(c_node)]] <- list(child = c_node, edge_indices = mask - 1L)
    }
    child_times <- node_times[unique_children + 1L]
    ord <- order(child_times, decreasing = TRUE)
    groups[as.character(unique_children[ord])]
}

.td_edges_by_child_then_parent_desc <- function(bp) {
    children <- bp$edge_children
    parents <- bp$edge_parents
    tc <- .td_tables(bp$ts)
    nt <- .td_nodes(tc)
    node_times <- nt$time
    unique_children <- unique(children)
    groups <- list()
    for (c_node in unique_children) {
        mask <- which(children == c_node)
        parent_times_for_child <- node_times[parents[mask] + 1L]
        edge_ord <- order(parent_times_for_child, decreasing = TRUE)
        groups[[as.character(c_node)]] <- list(child = c_node, edge_indices = mask[edge_ord] - 1L)
    }
    child_times <- node_times[unique_children + 1L]
    ord <- order(child_times, decreasing = TRUE)
    groups[as.character(unique_children[ord])]
}

.td_bp_inside_pass <- function(bp, standardize = TRUE, cache_inside = FALSE) {
    lik <- bp$lik
    inside <- .td_ntv_clone(bp$priors, grid_data = NaN, fixed_data = lik$identity_constant)
    num_edges <- bp$num_edges
    if (cache_inside) {
        g_i <- matrix(lik$identity_constant, nrow = num_edges, ncol = lik$grid_size)
    }
    denominator <- rep(NaN, .td_num_nodes(bp$ts))
    marginal_lik <- lik$identity_constant
    parent_groups <- .td_edges_by_parent_asc(bp)
    for (group in parent_groups) {
        parent <- group$parent
        edge_indices <- group$edge_indices
        if (parent %in% bp$fixed_set) 
            next
        val <- .td_ntv_get(bp$priors, parent)
        for (ei in edge_indices) {
            child <- bp$edge_children[ei + 1L]
            span <- bp$edge_spans[ei + 1L]
            spanfrac <- span/bp$spans[child + 1L]
            if (child %in% bp$fixed_set) {
                daughter_val <- .td_lik_scale_geometric(lik, spanfrac, .td_ntv_get(inside, child))
                edge_lik <- .td_lik_get_fixed(lik, daughter_val, ei)
            }
            else {
                inside_child <- .td_ntv_get(inside, child)
                if (all(is.nan(inside_child))) {
                  stop("The input tree sequence includes dangling nodes: ", "please simplify it")
                }
                daughter_val <- .td_lik_scale_geometric(lik, spanfrac, .td_lik_make_lower_tri(lik, inside_child))
                edge_lik <- .td_lik_get_inside(lik, daughter_val, ei)
            }
            val <- .td_lik_combine(lik, val, edge_lik)
            if (cache_inside) {
                g_i[ei + 1L, ] <- edge_lik
            }
        }
        denominator[parent + 1L] <- if (standardize) 
            max(val)
        else lik$identity_constant
        inside <- .td_ntv_set(inside, parent, .td_lik_ratio(lik, val, denominator[parent + 1L]))
        if (standardize) {
            marginal_lik <- .td_lik_combine(lik, marginal_lik, denominator[parent + 1L])
        }
    }
    for (key in names(bp$root_spans)) {
        root <- as.integer(key)
        span_when_root <- bp$root_spans[[key]]
        spanfrac <- span_when_root/bp$spans[root + 1L]
        root_val <- .td_lik_scale_geometric(lik, spanfrac, .td_ntv_get(inside, root))
        marginal_lik <- .td_lik_combine(lik, marginal_lik, .td_lik_marginalize(lik, root_val))
    }
    bp$inside <- inside
    bp$denominator <- denominator
    bp$marginal_lik <- marginal_lik
    if (cache_inside) {
        bp$g_i <- g_i/denominator[bp$edge_children + 1L]
        if (lik$log_space) {
            bp$g_i <- g_i - denominator[bp$edge_children + 1L]
        }
    }
    bp
}

.td_bp_outside_pass <- function(bp, standardize = FALSE, ignore_oldest_root = FALSE) {
    if (is.null(bp$inside)) {
        stop("You have not yet run the inside algorithm")
    }
    lik <- bp$lik
    outside <- .td_ntv_clone(bp$inside, grid_data = 0, probability_space = .td_LIN_GRID)
    for (key in names(bp$root_spans)) {
        root <- as.integer(key)
        span_when_root <- bp$root_spans[[key]]
        outside <- .td_ntv_set(outside, root, span_when_root/bp$spans[root + 1L])
    }
    outside <- .td_ntv_force_probability_space(outside, bp$inside$probability_space)
    num_nodes <- .td_num_nodes(bp$ts)
    child_groups <- .td_edges_by_child_desc(bp)
    for (group in child_groups) {
        child <- group$child
        edge_indices <- group$edge_indices
        if (child %in% bp$fixed_set) 
            next
        val <- rep(lik$identity_constant, lik$grid_size)
        for (ei in edge_indices) {
            parent <- bp$edge_parents[ei + 1L]
            if (ignore_oldest_root && parent == num_nodes - 1L) 
                next
            if (parent %in% bp$fixed_set) {
                stop("Fixed nodes cannot currently be parents in the TS")
            }
            span <- bp$edge_spans[ei + 1L]
            spanfrac <- span/bp$spans[child + 1L]
            if (!is.null(bp$g_i)) {
                inside_div_gi <- .td_lik_ratio(lik, .td_ntv_get(bp$inside, parent), bp$g_i[ei + 1L, ], 
                  div_0_null = TRUE)
            }
            else {
                daughter_val <- .td_lik_scale_geometric(lik, spanfrac, .td_lik_make_lower_tri(lik, .td_ntv_get(bp$inside, 
                  child)))
                edge_lik <- .td_lik_get_inside(lik, daughter_val, ei)
                cur_g_i <- .td_lik_ratio(lik, edge_lik, bp$denominator[child + 1L])
                inside_div_gi <- .td_lik_ratio(lik, .td_ntv_get(bp$inside, parent), cur_g_i, div_0_null = TRUE)
            }
            parent_val <- .td_lik_scale_geometric(lik, spanfrac, .td_lik_make_upper_tri(lik, .td_lik_combine(lik, 
                .td_ntv_get(outside, parent), inside_div_gi)))
            if (standardize) {
                parent_val <- .td_lik_ratio(lik, parent_val, max(parent_val))
            }
            edge_lik <- .td_lik_get_outside(lik, parent_val, ei)
            val <- .td_lik_combine(lik, val, edge_lik)
        }
        outside <- .td_ntv_set(outside, child, .td_lik_ratio(lik, val, bp$denominator[child + 1L]))
        if (standardize) {
            outside <- .td_ntv_set(outside, child, .td_lik_ratio(lik, val, max(val)))
        }
    }
    bp$outside <- outside
    bp$posterior_grid <- .td_ntv_clone(outside, grid_data = .td_lik_combine(lik, bp$inside$grid_data, 
        outside$grid_data), fixed_data = NaN)
    bp
}

.td_bp_outside_maximization <- function(bp, eps = 1e-06) {
    if (is.null(bp$inside)) {
        stop("You have not yet run the inside algorithm")
    }
    lik <- bp$lik
    num_nodes <- .td_num_nodes(bp$ts)
    maximized_node_times <- integer(num_nodes)
    all_children <- unique(bp$edge_children)
    mrcas <- setdiff(seq(0L, num_nodes - 1L), all_children)
    for (i in mrcas) {
        if (!(i %in% bp$fixed_set)) {
            maximized_node_times[i + 1L] <- which.max(.td_ntv_get(bp$inside, i)) - 1L
        }
    }
    child_groups <- .td_edges_by_child_then_parent_desc(bp)
    for (group in child_groups) {
        child <- group$child
        edge_indices <- group$edge_indices
        if (child %in% bp$fixed_set) 
            next
        result <- NULL
        youngest_par_index <- NULL
        for (k in seq_along(edge_indices)) {
            ei <- edge_indices[k]
            parent <- bp$edge_parents[ei + 1L]
            par_time_idx <- maximized_node_times[parent + 1L]
            parent_time <- lik$timepoints[par_time_idx + 1L]
            muts <- lik$mut_edges[ei + 1L]
            span <- lik$edge_spans[ei + 1L]
            if (k == 1L) {
                youngest_par_index <- par_time_idx
                dt <- parent_time - lik$timepoints[seq_len(youngest_par_index + 1L)] + eps
                lambda <- dt * lik$mut_rate * span
                if (lik$log_space) {
                  ll_mut <- dpois(muts, lambda, log = TRUE)
                  result <- ll_mut - max(ll_mut)
                }
                else {
                  ll_mut <- dpois(muts, lambda)
                  mx <- max(ll_mut)
                  result <- if (mx > 0) 
                    ll_mut/mx
                  else ll_mut
                }
            }
            else {
                cur_parent_index <- par_time_idx
                if (cur_parent_index < youngest_par_index) {
                  youngest_par_index <- cur_parent_index
                }
                dt <- parent_time - lik$timepoints[seq_len(youngest_par_index + 1L)] + eps
                lambda <- dt * lik$mut_rate * span
                if (lik$log_space) {
                  ll_mut <- dpois(muts, lambda, log = TRUE)
                  ll_norm <- ll_mut - max(ll_mut)
                  result[seq_len(youngest_par_index + 1L)] <- .td_lik_combine(lik, ll_norm, result[seq_len(youngest_par_index + 
                    1L)])
                }
                else {
                  ll_mut <- dpois(muts, lambda)
                  mx <- max(ll_mut)
                  ll_norm <- if (mx > 0) 
                    ll_mut/mx
                  else ll_mut
                  result[seq_len(youngest_par_index + 1L)] <- .td_lik_combine(lik, ll_norm, result[seq_len(youngest_par_index + 
                    1L)])
                }
            }
        }
        inside_val <- .td_ntv_get(bp$inside, child)[seq_len(youngest_par_index + 1L)]
        combined <- .td_lik_combine(lik, result[seq_len(youngest_par_index + 1L)], inside_val)
        maximized_node_times[child + 1L] <- which.max(combined) - 1L
    }
    bp$posterior_mean <- lik$timepoints[maximized_node_times + 1L]
    bp
}

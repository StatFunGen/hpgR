# Classroom R-list accessors for the tsdate routines from Daniel Nachun's tskitR.
# These adapters retain zero-based IDs and half-open genomic intervals.

.td_is <- function(object, class) {
    if (class == "TreeSequence") inherits(object, "hpgR_ancestry") else inherits(object, class)
}
.td_tables <- function(ts) ts
.td_nodes <- function(tc) tc$nodes
.td_edges <- function(tc) tc$edges
.td_num_nodes <- function(ts) nrow(ts$nodes)
.td_num_rows <- function(table) nrow(table)
.td_samples <- function(ts) ts$sampleIds
.td_tree <- function(ts) {
    x <- new.env(parent = emptyenv())
    x$ts <- ts
    x$cuts <- sort(unique(c(0, ts$sequenceLength, ts$edges$left, ts$edges$right)))
    x$index <- 1L
    x
}
.td_first <- function(tr) { tr$index <- 1L; invisible(TRUE) }
.td_next_tree <- function(tr) {
    if (tr$index >= length(tr$cuts) - 1L) return(FALSE)
    tr$index <- tr$index + 1L
    TRUE
}
.td_interval <- function(tr) tr$cuts[c(tr$index, tr$index + 1L)]
.td_num_trees <- function(ts) length(.td_tree(ts)$cuts) - 1L
.td_parent_node <- function(tr, node) {
    p <- .ancestryParent(tr$ts, mean(.td_interval(tr)))
    unname(p[as.character(node)])
}
.td_preorder <- function(tr) {
    p <- .ancestryParent(tr$ts, mean(.td_interval(tr)))
    unique(unlist(lapply(tr$ts$sampleIds, .ancestryPath, parent = p)))
}
.td_roots <- function(tr) {
    present <- .td_preorder(tr)
    p <- .ancestryParent(tr$ts, mean(.td_interval(tr)))
    present[p[as.character(present)] == -1L]
}
.td_num_samples_node <- function(tr, node) {
    p <- .ancestryParent(tr$ts, mean(.td_interval(tr)))
    sum(vapply(tr$ts$sampleIds, function(s) node %in% .ancestryPath(s, p), logical(1)))
}
.td_count_edge_mutations <- function(ts) {
    # Dan's .count_edge_mutations / tsdate_count_mutations: assign each
    # mutation to the incoming edge of its node at its genomic position.
    result <- integer(nrow(ts$edges))
    if (nrow(ts$mutations)) for (j in seq_len(nrow(ts$mutations))) {
        pos <- ts$sites$position[match(ts$mutations$site[j], ts$sites$id)]
        edge <- which(ts$edges$child == ts$mutations$node[j] &
                      ts$edges$left <= pos & pos < ts$edges$right)
        if (length(edge) == 1L) result[edge] <- result[edge] + 1L
    }
    result
}
.td_population_size_history <- function(population_size) {
    structure(list(size = .positiveScalar(population_size, "populationSize")),
              class = "PopulationSizeHistory")
}
.td_to_natural_timescale <- function(object, time) time * 2 * object$size
.td_to_coalescent_timescale <- function(object, time) time / (2 * object$size)

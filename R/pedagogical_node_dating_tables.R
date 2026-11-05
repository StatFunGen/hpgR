# Classroom R-list accessors for the tsdate routines from Daniel Nachun's tskitR.
# These adapters retain zero-based IDs and half-open genomic intervals.

tdIs <- function(object, class) {
    if (class == "TreeSequence") inherits(object, "hpgR_ancestry") else inherits(object, class)
}
tdTables <- function(ts) ts
tdNodes <- function(tc) tc$nodes
tdEdges <- function(tc) tc$edges
tdNumNodes <- function(ts) nrow(ts$nodes)
tdNumRows <- function(table) nrow(table)
tdSamples <- function(ts) ts$sampleIds
tdTree <- function(ts) {
    x <- new.env(parent = emptyenv())
    x$ts <- ts
    x$cuts <- sort(unique(c(0, ts$sequenceLength, ts$edges$left, ts$edges$right)))
    x$index <- 1L
    x
}
tdFirst <- function(tr) { tr$index <- 1L; invisible(TRUE) }
tdNextTree <- function(tr) {
    if (tr$index >= length(tr$cuts) - 1L) return(FALSE)
    tr$index <- tr$index + 1L
    TRUE
}
tdInterval <- function(tr) tr$cuts[c(tr$index, tr$index + 1L)]
tdNumTrees <- function(ts) length(tdTree(ts)$cuts) - 1L
tdParentNode <- function(tr, node) {
    p <- ancestryParent(tr$ts, mean(tdInterval(tr)))
    unname(p[as.character(node)])
}
tdPreorder <- function(tr) {
    p <- ancestryParent(tr$ts, mean(tdInterval(tr)))
    unique(unlist(lapply(tr$ts$sampleIds, ancestryPath, parent = p)))
}
tdRoots <- function(tr) {
    present <- tdPreorder(tr)
    p <- ancestryParent(tr$ts, mean(tdInterval(tr)))
    present[p[as.character(present)] == -1L]
}
tdNumSamplesNode <- function(tr, node) {
    p <- ancestryParent(tr$ts, mean(tdInterval(tr)))
    sum(vapply(tr$ts$sampleIds, function(s) node %in% ancestryPath(s, p), logical(1)))
}
tdCountEdgeMutations <- function(ts) {
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
tdPopulationSizeHistory <- function(population_size) {
    structure(list(size = positiveScalar(population_size, "populationSize")),
              class = "PopulationSizeHistory")
}
tdToNaturalTimescale <- function(object, time) time * 2 * object$size
tdToCoalescentTimescale <- function(object, time) time / (2 * object$size)

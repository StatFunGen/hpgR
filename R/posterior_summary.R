# From Daniel Nachun's demographyTools: .posterior_summary.

#' Summarize samples from a posterior distribution
#' @param samples Numeric matrix, draws in rows and parameters in columns.
#' @param transform_fn Optional transformation applied to each row.
#' @param probs Quantile probabilities; first and last define the reported interval.
#' @return Data frame of posterior means, medians, standard deviations and quantiles.
#' @export
posteriorSummary <- function(samples, transform_fn = NULL,
                             probs = c(0.025, 0.5, 0.975)) {
    if (!is.null(transform_fn)) {
        samples <- t(apply(samples, 1, transform_fn))
    }
    n_params <- ncol(samples)
    data.frame(
        mean = colMeans(samples),
        median = apply(samples, 2, stats::median),
        sd = apply(samples, 2, stats::sd),
        lower = apply(samples, 2, stats::quantile, probs = probs[1]),
        upper = apply(samples, 2, stats::quantile, probs = probs[length(probs)])
    )
}

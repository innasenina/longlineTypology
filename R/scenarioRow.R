#' Record one scenario's clustering result as a single row
#'
#' Captures the two things the manuscript never reports -- the number of
#' retained principal components, and the distribution of optimal K across
#' replicates -- alongside the metrics it does report, so that scenarios can be
#' compared with the confounds visible.
#'
#' @param scenario   scenario string, e.g. "yba_lat".
#' @param pca_res    the object returned by customPCA().
#' @param freq_table table() of optimal K across the replicates.
#' @param n_sub      subsample size used per replicate.
#' @param silhouette,nrmse optional, if computed elsewhere.
#' @return a one-row data frame; rbind() successive calls together.
#' @export
scenarioRow <- function(scenario, pca_res, freq_table, n_sub = NA,
						silhouette = NA, nrmse = NA) {
	ev  <- pca_res$sdev^2
	pcs <- pca_res$no_var
	data.frame(
		scenario    = scenario,
		n_sub       = n_sub,
		n_vars      = length(ev),
		PCs         = pcs,
		var_expl    = round(sum(ev[seq_len(pcs)]) / sum(ev), 3),
		K           = as.integer(names(freq_table)[which.max(freq_table)]),
		K_stable    = max(freq_table) == sum(freq_table),
		K_dist      = paste(names(freq_table), as.integer(freq_table),
							sep = "x", collapse = " "),
		eigen_rel   = paste(round(ev / ev[1], 4), collapse = " "),
		silhouette  = silhouette,
		nrmse       = nrmse,
		stringsAsFactors = FALSE)
}

#' Append a scenario row to a running results table (creates it if absent)
#' @export
addScenario <- function(results, ...) {
	row <- scenarioRow(...)
	if (is.null(results) || !nrow(results)) row else rbind(results, row)
}

#' One-call profile of a clustering, with clusters labelled by dominant species
#'
#' Replaces the five ad-hoc tapply() calls, and labels each cluster by which
#' fraction is largest so downstream code never hard-codes a cluster number.
#'
#' @param df    the data frame. Fraction columns are required; `mean_hbf` and
#'   `mean_len` are profiled if present.
#'
#'   Note on `len`: it is only as good as the length fill behind it. Under MICE
#'   the values are imputed FROM species composition, so a length difference
#'   between clusters partly restates the composition difference that defined
#'   them -- see B2. `len_n` gives the non-missing count per cluster so the
#'   coverage is visible. The paper reports mean length per fishery in Table 4
#'   (117.7 / 118.9 / 113.3 cm), so this column is directly comparable to it.
#' @param cl    cluster vector, same length as nrow(df).
#' @param share optional vector, e.g. yba_share_w.
#' @param zero  optional logical/0-1 vector, e.g. total_n == 0.
#' @param vars  fraction columns to profile.
#' @return a data frame, one row per cluster, with a `dominant` label.
#' @export
clusterProfile <- function(df, cl, share = NULL, zero = NULL,
						   vars = c("yft_fraction", "bet_fraction",
									"alb_fraction", "oth_fraction")) {
	cl <- as.factor(cl)
	m  <- sapply(split(df[, vars, drop = FALSE], cl),
				 function(d) sapply(d, stats::median, na.rm = TRUE))
	out <- data.frame(cluster = levels(cl), t(round(m, 3)),
					  stringsAsFactors = FALSE, row.names = NULL)
	out$dominant <- sub("_fraction", "", rownames(m)[apply(m, 2L, which.max)])
	if ("mean_hbf" %in% names(df))
		out$hbf <- round(tapply(df$mean_hbf, cl, stats::median, na.rm = TRUE), 1)
		out$hbf_mean <- round(tapply(df$mean_hbf, cl, mean, na.rm = TRUE), 1)
	if ("mean_len" %in% names(df)) {
		out$len   <- round(tapply(df$mean_len, cl, stats::median, na.rm = TRUE), 1)
		out$len_n <- as.integer(tapply(!is.na(df$mean_len), cl, sum))
	}
	if (!is.null(share))
		out$target_share <- round(tapply(share, cl, stats::median, na.rm = TRUE), 3)
	if (!is.null(zero))
		out$zero_rate <- round(tapply(as.numeric(zero), cl, mean), 4)
	out$size <- round(as.numeric(prop.table(table(cl))), 3)
	out
}

#' Which cluster(s) are the non-target fishery?
#'
#' Returns cluster labels, not positions -- use with `cl %in% nonTarget(prof)`.
#' A cluster qualifies if `oth` is its dominant fraction, or if its target share
#' is below `share_max`.
#'
#' **Results are ordered by target share, purest first**, because more than one
#' cluster commonly qualifies and they are not interchangeable: at K = 5 the
#' `sp` scenario yields a *pure* non-target group (target share 0.15, droppable
#' from the model) and a *mixed* one (0.67, a real fishery that lands bigeye).
#' Use `nonTarget(prof)[1]` for the pure group; use the whole vector only when
#' you mean their union, which reproduces the coarser K = 4 cluster.
#'
#' @param prof output of clusterProfile().
#' @param share_max clusters below this target share qualify regardless of
#'   which fraction dominates.
#' @param purest if TRUE, return only the single lowest-target-share cluster.
#' @export
nonTarget <- function(prof, share_max = 0.6, purest = FALSE) {
	i <- prof$dominant == "oth"
	if (!is.null(prof$target_share)) i <- i | prof$target_share < share_max
	out <- prof[i, , drop = FALSE]
	if (!is.null(out$target_share)) out <- out[order(out$target_share), ]
	if (purest) out$cluster[1L] else out$cluster
}

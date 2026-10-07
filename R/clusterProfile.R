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
#' @param date_col name of the date column used for the season columns; by
#'   default the first of "date" or "ymd" present. If none is found, the season
#'   columns are skipped.
#' @param effort_col name of the effort column for the effort-weighted season
#'   columns; skipped if absent.
#' @param lat_col name of the latitude column that decides the hemisphere.
#' @return a data frame, one row per cluster, with a `dominant` label. With
#'   date and latitude columns it also has `winter`, `spring`, `summer`,
#'   `autumn`: the share of the cluster's records in each LOCAL season (see
#'   [seasonOf()]), summing to 1 across a row; and, with an effort column,
#'   `winter_e`..`autumn_e`: the share of the cluster's effort in each season.
#'   About 0.25 each means the fishery operates year-round; compare against
#'   [seasonOverall()] if the data themselves are uneven by season. Records
#'   answer "is the fishery present", effort "how much it fishes".
#' @export
clusterProfile <- function(df, cl, share = NULL, zero = NULL,
						   vars = c("yft_fraction", "bet_fraction",
									"alb_fraction", "oth_fraction"),
						   date_col = NULL, effort_col = "E",
						   lat_col = "latitude") {
	cl <- as.factor(cl)
	m  <- sapply(split(df[, vars, drop = FALSE], cl),
				 function(d) sapply(d, stats::median, na.rm = TRUE))
	out <- data.frame(cluster = levels(cl), t(round(m, 3)),
					  stringsAsFactors = FALSE, row.names = NULL)
	out$dominant <- sub("_fraction", "", rownames(m)[apply(m, 2L, which.max)])
	if ("mean_hbf" %in% names(df)) {
		out$hbf <- round(tapply(df$mean_hbf, cl, stats::median, na.rm = TRUE), 1)
		out$hbf_mean <- round(tapply(df$mean_hbf, cl, mean, na.rm = TRUE), 1)
	}
	if ("mean_len" %in% names(df)) {
		out$len   <- round(tapply(df$mean_len, cl, stats::median, na.rm = TRUE), 1)
		out$len_n <- as.integer(tapply(!is.na(df$mean_len), cl, sum))
	}
	if (!is.null(share))
		out$target_share <- round(tapply(share, cl, stats::median, na.rm = TRUE), 3)
	if (!is.null(zero))
		out$zero_rate <- round(tapply(as.numeric(zero), cl, mean), 4)
	q <- seasonOf(df, date_col, lat_col)
	if (!is.null(q)) {
		qs <- round(prop.table(table(cl, q), 1L), 3)
		for (j in levels(q)) out[[j]] <- as.numeric(qs[, j])
		if (effort_col %in% names(df)) {
			qe <- round(effortByQuarter(df[[effort_col]], cl, q), 3)
			for (j in levels(q)) out[[paste0(j, "_e")]] <- as.numeric(qe[, j])
		}
	}
	out$size <- round(as.numeric(prop.table(table(cl))), 3)
	out
}

## Calendar quarter (factor 1..4) from the first date column found, or NULL.
quarterOf <- function(df, date_col = NULL) {
	if (is.null(date_col)) date_col <- intersect(c("date", "ymd"), names(df))[1L]
	if (is.na(date_col) || !date_col %in% names(df)) return(NULL)
	m <- as.integer(format(as.Date(df[[date_col]]), "%m"))
	factor((m - 1L) %/% 3L + 1L, levels = 1:4)
}

#' Local (hemisphere-relative) meteorological season of each record
#'
#' Winter = Dec-Feb north of the equator and Jun-Aug south of it; spring =
#' Mar-May / Sep-Nov; summer = Jun-Aug / Dec-Feb; autumn = Sep-Nov / Mar-May.
#' Latitude 0 counts as northern.
#' @param df data frame with a date column and `lat_col`.
#' @param date_col as in [clusterProfile()].
#' @param lat_col latitude column.
#' @return factor with levels winter, spring, summer, autumn; NULL if either
#'   column is missing.
#' @export
seasonOf <- function(df, date_col = NULL, lat_col = "latitude") {
	if (is.null(date_col)) date_col <- intersect(c("date", "ymd"), names(df))[1L]
	if (is.na(date_col) || !date_col %in% names(df) || !lat_col %in% names(df))
		return(NULL)
	m <- as.integer(format(as.Date(df[[date_col]]), "%m"))
	s <- (m %% 12L) %/% 3L                   # Dec-Feb 0, Mar-May 1, Jun-Aug 2, Sep-Nov 3
	s <- ifelse(df[[lat_col]] < 0, (s + 2L) %% 4L, s)   # southern: shift by six months
	factor(c("winter", "spring", "summer", "autumn")[s + 1L],
		   levels = c("winter", "spring", "summer", "autumn"))
}

## Share of each group's effort by quarter or season (rows sum to 1). Records with
## missing effort are left out of both numerator and denominator.
effortByQuarter <- function(e, g, q) {
	ok <- !is.na(e)
	s  <- tapply(e[ok], list(g[ok], q[ok]), sum)
	s[is.na(s)] <- 0
	s / rowSums(s)
}

#' Seasonal shares over all records, the baseline for the season profile columns
#'
#' @param df the data frame given to [clusterProfile()].
#' @param date_col,effort_col,lat_col as in [clusterProfile()].
#' @return a matrix with columns winter..autumn and a `records` row, plus an
#'   `effort` row when the effort column is present; NULL without date or
#'   latitude.
#' @export
seasonOverall <- function(df, date_col = NULL, effort_col = "E",
						  lat_col = "latitude") {
	q <- seasonOf(df, date_col, lat_col)
	if (is.null(q)) return(NULL)
	out <- rbind(records = as.numeric(prop.table(table(q))))
	if (effort_col %in% names(df))
		out <- rbind(out, effort = effortByQuarter(df[[effort_col]],
												   rep(1L, nrow(df)), q)[1L, ])
	colnames(out) <- levels(q)
	round(out, 3)
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

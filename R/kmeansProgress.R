#' k-means over nstart restarts, with a progress bar and a local-optimum report
#'
#' [stats::kmeans()] runs its `nstart` restarts inside compiled code, so nothing
#' can report progress from within it. This drives the restarts from R instead:
#' `nstart` calls with `nstart = 1`, keeping the one with the lowest
#' `tot.withinss`. That is exactly what `kmeans()` does internally, so the
#' procedure is unchanged.
#'
#' **It is not bit-identical to `kmeans(nstart = N)`**, because the restarts are
#' drawn from the RNG in a different order. Same procedure, same distribution,
#' different draw -- so a run made this way will not reproduce a stored result
#' made the other way. Set `progress = FALSE` to take the single-call path and
#' reproduce older runs exactly.
#'
#' The report afterwards is the part worth having. With 20,000+ points and
#' K = 5 the objective has competing optima a couple of percent apart, and most
#' seeds land on the worse one, so knowing how many restarts reached the best
#' solution says whether `nstart` is high enough.
#'
#' @param x,centers,iter.max,algorithm as for [stats::kmeans()].
#' @param nstart number of restarts.
#' @param progress show the bar and the report. FALSE reverts to one
#'   `kmeans(nstart = nstart)` call.
#' @param tol relative gap within which a restart counts as having reached the
#'   best solution.
#' @param material relative spread below which the restarts are treated as the
#'   same solution, and the count of restarts reaching the best is suppressed as
#'   uninformative.
#' @return the best [stats::kmeans()] object, with `tot.withinss` from every
#'   restart attached as attribute `"starts"`.
#' @family clustering utilities
#' @export
kmeansProgress <- function(x, centers, iter.max = 10, nstart = 1,
						   algorithm = "Hartigan-Wong",
						   progress = TRUE, tol = 1e-6, material = 1e-3) {
	if (!progress || nstart <= 1)
		return(stats::kmeans(x, centers = centers, iter.max = iter.max,
							 nstart = nstart, algorithm = algorithm))

	pb <- utils::txtProgressBar(min = 0, max = nstart, style = 3)
	on.exit(close(pb), add = TRUE)

	best <- NULL; w <- rep(NA_real_, nstart)
	for (i in seq_len(nstart)) {
		fit <- try(stats::kmeans(x, centers = centers, iter.max = iter.max,
								 nstart = 1L, algorithm = algorithm),
				   silent = TRUE)
		if (!inherits(fit, "try-error")) {
			w[i] <- fit$tot.withinss
			if (is.null(best) || fit$tot.withinss < best$tot.withinss) best <- fit
		}
		utils::setTxtProgressBar(pb, i)
	}
	close(pb); on.exit()
	if (is.null(best)) stop("every restart failed")

	ok     <- stats::na.omit(w)
	spread <- max(ok) / min(ok) - 1
	hit    <- sum(ok <= min(ok) * (1 + tol))
	## Report the SPREAD first. The count of restarts reaching the best is only
	## meaningful when the solutions actually differ: at K = 3 in two dimensions
	## the spread is ~1e-4, so "2/30 reached the best" is counting rounding, not
	## competing optima, and reads far more alarming than it is.
	cat(sprintf("\n%d restarts, spread %+.3f%% (tot.withinss %.6g)",
				length(ok), 100 * spread, min(ok)))
	if (length(ok) < nstart) cat(sprintf(", %d failed", nstart - length(ok)))
	if (spread < material) {
		cat(" -- no competing optima, nstart is not binding here\n")
	} else {
		cat(sprintf("\n  %d/%d restarts reached the best solution\n", hit, nstart))
		if (hit == 1L)
			cat("  only one restart found it -- raise nstart before trusting this partition\n")
	}
	utils::flush.console()

	attr(best, "starts") <- w
	best
}

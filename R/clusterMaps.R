#' Grid catch, effort and CPUE per cluster onto a lon/lat mesh
#'
#' Grids catch and effort separately onto a regular lon/lat mesh for each
#' cluster and returns CPUE as the *ratio of sums*, not the mean of per-record
#' ratios. That matters: a cell holding one lucky set with tiny effort would
#' otherwise dominate its own average. This is the same convention
#' [summaryClusters()] uses, where CPUE is a per-(cluster, day) sum of catch
#' over a sum of effort.
#'
#' Cells with little effort are unstable and are masked by `min_effort` rather
#' than plotted as noise. The default keeps the map honest at 5 degrees; lower
#' it if you grid finer.
#'
#' @param df   data frame with `longitude`, `latitude`, the effort column, and
#'   the species count columns.
#' @param cl   cluster vector, same length as `nrow(df)`.
#' @param species one of `"yft"`, `"bet"`, `"alb"`, `"oth"` (oth + skj),
#'   `"yba"` (the three tunas) or `"all"`.
#' @param res  grid size in degrees. 1 matches the data but is noisy over the
#'   whole basin; 5 is readable.
#' @param effort_col name of the effort column, hundreds of hooks.
#' @param min_effort cells with less total effort than this are set to NA.
#' @param xlim,ylim longitude (0-360 convention) and latitude range; taken from
#'   the data when NULL.
#'
#' @return list with `cpue` (one matrix per cluster), `catch`, `effort`, `hbf`
#'   (effort-weighted mean hooks between floats, `NULL` if the column is absent),
#'   `effort_raw` (the same effort grid **before** the `min_effort` mask -- shares
#'   must be computed from this, see [plotClusterDominance()]), the
#'   `lon`/`lat` cell centres, and the settings used. Matrices are
#'   `length(lon)` by `length(lat)`, ready for [image()].
#' @family clustering utilities
#' @export
clusterGrid <- function(df, cl, species = "yft", res = 5, effort_col = "E",
					 min_effort = 100, xlim = NULL, ylim = NULL) {
	stopifnot(nrow(df) == length(cl))
	need <- c("longitude", "latitude", effort_col)
	miss <- need[!need %in% names(df)]
	if (length(miss)) stop("missing column(s): ", paste(miss, collapse = ", "))

	catch <- switch(species,
		yft = df$yft_n, bet = df$bet_n, alb = df$alb_n,
		oth = df$oth_n + df$skj_n,
		yba = ybaN(df), all = allN(df),
		stop("species must be one of yft, bet, alb, oth, yba, all"))
	eff <- df[[effort_col]]

	## 0-360 so the Pacific is contiguous
	lon <- df$longitude %% 360
	lat <- df$latitude
	ok  <- is.finite(lon) & is.finite(lat) & is.finite(catch) & is.finite(eff) & eff > 0
	if (!any(ok)) stop("no usable records")

	if (is.null(xlim)) xlim <- range(lon[ok])
	if (is.null(ylim)) ylim <- range(lat[ok])
	xb <- seq(floor(xlim[1] / res) * res, ceiling(xlim[2] / res) * res, by = res)
	yb <- seq(floor(ylim[1] / res) * res, ceiling(ylim[2] / res) * res, by = res)
	nx <- length(xb) - 1L; ny <- length(yb) - 1L

	ix <- findInterval(lon, xb, rightmost.closed = TRUE)
	iy <- findInterval(lat, yb, rightmost.closed = TRUE)
	ok <- ok & ix >= 1L & ix <= nx & iy >= 1L & iy <= ny

	clf  <- factor(cl[ok])
	## factor() drops unused levels, so a cluster whose records were all removed
	## by the filters above would vanish without trace. Say so.
	if (is.factor(cl) && nlevels(clf) < nlevels(cl))
		warning("cluster(s) dropped by filtering and absent from the grid: ",
				paste(setdiff(levels(cl), levels(clf)), collapse = ", "))
	cell <- (iy[ok] - 1L) * nx + ix[ok]          # column-major, matches matrix()
	lev  <- levels(clf)

	## rowsum() over a (cluster, cell) key does both grids in one pass
	## HBF is accumulated as sum(hbf * E) and sum(E over records WITH an hbf), so
	## the cell value is the effort-weighted mean -- the depth most hooks were
	## fished at, not the average over records. Records without a measured hbf
	## are excluded from both sums, otherwise their effort would drag the mean
	## down while contributing no depth.
	has_h <- "mean_hbf" %in% names(df)
	hv <- if (has_h) df$mean_hbf[ok] else rep(NA_real_, sum(ok))
	fin <- is.finite(hv)
	key <- paste(as.integer(clf), cell, sep = "\r")
	m   <- rowsum(cbind(catch  = catch[ok], effort = eff[ok],
						hbf_we = ifelse(fin, hv * eff[ok], 0),
						hbf_e  = ifelse(fin, eff[ok], 0)),
				  group = key, reorder = FALSE)
	kk  <- strsplit(rownames(m), "\r", fixed = TRUE)
	ki  <- as.integer(vapply(kk, `[`, "", 1L))
	kc  <- as.integer(vapply(kk, `[`, "", 2L))

	blank <- function() matrix(NA_real_, nx, ny)
	out_c <- out_e <- out_p <- out_h <- out_r <-
		stats::setNames(vector("list", length(lev)), lev)
	for (j in seq_along(lev)) {
		s <- ki == j
		C <- blank(); E <- blank(); HW <- blank(); HE <- blank()
		C[kc[s]]  <- m[s, "catch"];  E[kc[s]]  <- m[s, "effort"]
		HW[kc[s]] <- m[s, "hbf_we"]; HE[kc[s]] <- m[s, "hbf_e"]
		out_r[[j]] <- E                      # unmasked, for dominance shares
		E[!is.na(E) & E < min_effort] <- NA_real_
		HE[!is.na(HE) & HE < min_effort] <- NA_real_
		HE[!is.na(HE) & HE == 0] <- NA_real_
		out_c[[j]] <- C; out_e[[j]] <- E; out_p[[j]] <- C / E
		out_h[[j]] <- if (has_h) HW / HE else NULL
	}
	list(cpue = out_p, catch = out_c, effort = out_e, effort_raw = out_r,
		 hbf = if (has_h) out_h else NULL,
		 lon = xb[-length(xb)] + res / 2, lat = yb[-length(yb)] + res / 2,
		 breaks_lon = xb, breaks_lat = yb,
		 species = species, res = res, min_effort = min_effort)
}

#' Draw the maps returned by clusterGrid()
#'
#' The field drawn is chosen by `what`, so check the colour-bar label: `"effort"`
#' and `"cpue"` are very different maps of the same fishery.
#'
#' One panel per cluster on a **shared** colour scale, which is the whole point:
#' panels scaled individually look alike no matter how different they are. CPUE
#' is heavily right-skewed, so the default scale is log.
#'
#' Coastlines are drawn if the `maps` package happens to be installed, and
#' silently skipped otherwise, so this adds no dependency.
#'
#' @param g output of [clusterGrid()].
#' @param what which field to draw. `"cpue"` is where a fishery caught well,
#'   `"effort"` is where it actually operated -- use `"effort"` to judge spatial
#'   footprint, since a cluster with little effort in a region still shows a CPUE
#'   there. `"share"` gives each cluster's fraction of the cell's total effort,
#'   which is the direct picture of spatial mixing between clusters.
#' @param log log10 colour scale. Default TRUE, except for `"share"`.
#' @param zlim shared colour limits; quantiles of the pooled data when NULL.
#' @param probs quantiles used for `zlim`, to keep outliers from eating the ramp.
#' @param labels optional names for the panels, **in cluster order**. Cluster
#'   indices are arbitrary and change between runs, so prefer `prof=`.
#' @param prof a [clusterProfile()] table; labels are then derived from it.
#' @param mfrow panel layout; chosen automatically when NULL.
#' @param col colour ramp.
#' @param asp panel aspect ratio. NA (default) fills the panel with exactly the
#'   data extent. `asp = 1` gives true geometry but base graphics achieves it by
#'   *expanding the window*, not shrinking the panel, so the map is then padded
#'   with neighbouring ocean and coastline.
#' @family clustering utilities
#' @export
plotClusterMaps <- function(g, what = c("cpue", "effort", "catch", "share", "hbf"),
						 log = NULL, zlim = NULL, probs = c(0.02, 0.98),
						 labels = NULL, prof = NULL, mfrow = NULL, col = NULL, asp = NA) {
	what <- match.arg(what)
	z <- switch(what, cpue = g$cpue, effort = g$effort, catch = g$catch,
				hbf = g$hbf,
				share = local({
					tot <- Reduce(function(a, b) {
						a[is.na(a)] <- 0; b[is.na(b)] <- 0; a + b }, g$effort)
					tot[tot == 0] <- NA_real_
					lapply(g$effort, function(m) m / tot)
				}))
	if (what == "hbf" && is.null(z))
		stop("no hbf grid: clusterGrid() found no `mean_hbf` column in the data")
	## hbf spans roughly 2-50 and is read directly as a gear setting, so it is
	## drawn on a linear scale; the others are heavily right-skewed.
	if (is.null(log))  log  <- !(what %in% c("share", "hbf"))
	if (what == "share" && is.null(zlim)) zlim <- c(0, 1)
	if (is.null(col))
		col <- grDevices::hcl.colors(64, "YlGnBu", rev = TRUE)
	tr <- function(x) if (log) log10(x) else x

	pool <- unlist(lapply(z, function(m) m[is.finite(m) & m > 0]), use.names = FALSE)
	if (!length(pool)) stop("no finite positive CPUE to plot")
	if (is.null(zlim)) zlim <- stats::quantile(tr(pool), probs, na.rm = TRUE)
	else zlim <- tr(zlim)

	n <- length(z)
	if (is.null(mfrow)) { nc <- ceiling(sqrt(n + 1)); mfrow <- c(ceiling((n + 1) / nc), nc) }
	if (!is.null(prof)) {
		if (nrow(prof) != length(z))
			stop("prof has ", nrow(prof), " rows but the grid holds ", length(z),
				 " clusters -- are they from the same run?")
		## Only a genuine re-ordering is a problem. Once the cluster vector has
		## been named -- factor(cl, labels = clusterLabels(prof)) -- prof$cluster
		## holds ids and the grid holds names; different representations of the
		## same thing, and nothing to warn about.
		pc <- as.character(prof$cluster); nz <- names(z)
		if (!is.null(prof$cluster) && !identical(pc, nz) && setequal(pc, nz))
			warning("prof$cluster (", paste(pc, collapse = ","),
					") is in a different order from the grid's cluster levels (",
					paste(nz, collapse = ","), ")")
	}
	if (is.null(labels) && !is.null(prof)) labels <- clusterLabels(prof)
	if (is.null(labels)) labels <- names(z)

	op <- graphics::par(mfrow = mfrow, mar = c(2.5, 2.5, 2, 1), mgp = c(1.5, 0.4, 0),
						tcl = -0.25, cex.axis = 0.8)
	on.exit(graphics::par(op), add = TRUE)
	has_maps <- requireNamespace("maps", quietly = TRUE)

	for (j in seq_len(n)) {
		m <- tr(z[[j]])
		m[is.finite(m) & m < zlim[1]] <- zlim[1]
		m[is.finite(m) & m > zlim[2]] <- zlim[2]
		graphics::image(g$lon, g$lat, m, col = col, zlim = zlim,
						xlab = "", ylab = "", main = labels[j], useRaster = TRUE, asp = asp)
		if (has_maps)
			maps::map("world2", add = TRUE, interior = FALSE,
					  col = "grey30", lwd = 0.5)
		graphics::box()
	}

	## colour bar in the spare panel: a thin strip, not a stretched image
	graphics::par(mar = c(3, 2, 3, 2))
	b <- seq(zlim[1], zlim[2], length.out = length(col) + 1)
	graphics::plot.new()
	graphics::plot.window(xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i")
	graphics::rect(0.10, b[-length(b)], 0.40, b[-1], col = col, border = NA)
	graphics::rect(0.10, zlim[1], 0.40, zlim[2], border = "grey30", lwd = 0.6)
	at <- pretty(b, 5); at <- at[at >= zlim[1] & at <= zlim[2]]
	graphics::axis(4, at = at, pos = 0.40, las = 1, cex.axis = 0.8, tcl = -0.2,
				   labels = if (log) formatC(10^at, format = "g", digits = 3) else at)
	ttl <- switch(what,
		cpue   = paste0(g$species, " CPUE (fish per 100 hooks)"),
		effort = "effort (hundred hooks)",
		catch  = paste0(g$species, " catch (fish)"),
		share  = "share of cell effort",
		hbf    = "hooks between floats")
	graphics::mtext(ttl, side = 3,
					line = 0.5, cex = 0.8, adj = 0)
	graphics::mtext(sprintf("%d\u00b0 cells, effort >= %g", g$res, g$min_effort),
					side = 1, line = 0.5, cex = 0.7, adj = 0, col = "grey30")
	invisible(zlim)
}

#' Which fishery owns each cell, and how cleanly
#'
#' Maps the cluster holding the largest share of each cell's effort. Spatially
#' distinct fisheries give clean blocks of colour; a fishery defined only by
#' catch composition can be spread through the same water as every other one,
#' which this shows at a glance and the CPUE panels do not.
#'
#' The printed `separation` is the effort-weighted mean of the winning share.
#' It runs from 1/K, meaning every cell is split evenly between all K clusters,
#' to 1, meaning every cell belongs to a single one.
#'
#' @param g output of [clusterGrid()].
#' @param labels legend names, **in cluster order**. Prefer `prof=`.
#' @param prof a [clusterProfile()] table; labels are then derived from it.
#' @param col qualitative palette, one colour per cluster.
#' @param min_share cells whose winner holds less than this are left blank. At 0
#'   (the default) a cell is coloured whether its winner holds 100% or a bare
#'   plurality, so the map shows who wins and not how clearly; `share` is returned
#'   so that can be mapped separately.
#' @param min_cell_effort cells with less TOTAL effort than this are excluded.
#'   Defaults to the grid's `min_effort`.
#' @return invisibly, a list with the `winner` and `share` matrices and the
#'   scalar `separation`.
#' @family clustering utilities
#' @export
plotClusterDominance <- function(g, labels = NULL, prof = NULL, col = NULL,
								 min_share = 0, min_cell_effort = NULL) {
	## Shares MUST come from unmasked effort. `min_effort` in clusterGrid() is a
	## per-fishery mask written for CPUE, where a ratio from tiny effort is
	## unstable. Applied here it deletes small fisheries from the denominator, so
	## a genuinely contested cell -- say 120/90/80 -- is recorded as exclusively
	## one fishery's, share 1.000 instead of 0.414. The bias grows with K, since
	## more fisheries means less effort each and more of them fall below the cut,
	## which would make a higher-K partition look more cleanly separated purely
	## as an artefact. The question dominance asks is about the CELL: was it
	## fished enough to say who owns it? So the threshold belongs on the total.
	raw <- if (!is.null(g$effort_raw)) g$effort_raw else g$effort
	E <- lapply(raw, function(m) { m[is.na(m)] <- 0; m })
	K <- length(E)
	tot <- Reduce(`+`, E)
	if (is.null(min_cell_effort)) min_cell_effort <- g$min_effort
	tot[tot < min_cell_effort] <- 0
	arr <- array(unlist(E), dim = c(dim(E[[1]]), K))
	win <- apply(arr, c(1, 2), which.max)
	shr <- apply(arr, c(1, 2), max) / tot
	win[!is.finite(shr) | tot == 0 | shr < min_share] <- NA
	shr[!is.finite(shr) | tot == 0] <- NA

	if (is.null(col))    col    <- grDevices::hcl.colors(K, "Dark 3")
	if (!is.null(prof)) {
		if (nrow(prof) != length(g$cpue))
			stop("prof has ", nrow(prof), " rows but the grid holds ", length(g$cpue),
				 " clusters -- are they from the same run?")
		## Only a genuine re-ordering is a problem. Once the cluster vector has
		## been named -- factor(cl, labels = clusterLabels(prof)) -- prof$cluster
		## holds ids and the grid holds names; different representations of the
		## same thing, and nothing to warn about.
		pc <- as.character(prof$cluster); nz <- names(g$cpue)
		if (!is.null(prof$cluster) && !identical(pc, nz) && setequal(pc, nz))
			warning("prof$cluster (", paste(pc, collapse = ","),
					") is in a different order from the grid's cluster levels (",
					paste(nz, collapse = ","), ")")
	}
	if (is.null(labels) && !is.null(prof)) labels <- clusterLabels(prof)
	if (is.null(labels)) labels <- names(g$cpue)
	op <- graphics::par(mar = c(2.5, 2.5, 2, 1), mgp = c(1.5, 0.4, 0), tcl = -0.25,
						cex.axis = 0.8)
	on.exit(graphics::par(op), add = TRUE)
	graphics::image(g$lon, g$lat, win, col = col, zlim = c(0.5, K + 0.5),
					xlab = "", ylab = "", useRaster = TRUE,
					main = "dominant fishery by cell effort")
	if (requireNamespace("maps", quietly = TRUE))
		maps::map("world2", add = TRUE, interior = FALSE, col = "grey30", lwd = 0.5)
	graphics::box()
	graphics::legend("bottomleft", legend = labels, fill = col, bty = "n",
					 cex = 0.7, bg = "white")

	sep <- stats::weighted.mean(shr, tot, na.rm = TRUE)
	graphics::mtext(sprintf(
		"separation %.2f  (even mixing = %.2f, exclusive = 1)   |   %g\u00b0 cells, cell effort >= %g, min_share = %g",
		sep, 1 / K, g$res, min_cell_effort, min_share),
		side = 1, line = 1.3, cex = 0.7, adj = 0)
	invisible(list(winner = win, share = shr, separation = sep))
}

#' How much of each variable survived the PCA truncation
#'
#' A variable can be named as a clustering variable and still play almost no
#' part, if its variance loads onto components the variance threshold discards.
#' On z-scored data each variable has unit variance, so the sum of its squared
#' loadings weighted by the eigenvalues, over the retained components, is the
#' fraction of it that reached the k-means step.
#'
#' @param pca_res the object returned by [customPCA()].
#' @param n_comp number of **PCA components** retained -- not a number of
#'   clusters. Defaults to `pca_res$no_var`, whatever the run actually used.
#'   Named to match [customPCA()]; `k` and `K` mean clusters everywhere else in
#'   this package.
#' @return named vector, fraction of each variable's variance retained.
#' @family clustering utilities
#' @export
varRetained <- function(pca_res, n_comp = pca_res$no_var) {
	L   <- pca_res$rotation
	ev  <- pca_res$sdev^2
	if (n_comp > ncol(L))
		stop("n_comp = ", n_comp, " exceeds the ", ncol(L), " components available")
	tot <- rowSums(sweep(L^2, 2, ev, `*`))                  # = 1 if z-scored
	kept<- rowSums(sweep(L[, seq_len(n_comp), drop = FALSE]^2, 2,
					     ev[seq_len(n_comp)], `*`))
	round(kept / tot, 3)
}

#' Fishery names in cluster order, read off a profile
#'
#' Cluster indices are arbitrary: k-means numbers them by whatever the starting
#' centres happened to be, so a hand-written `labels` vector silently rotates
#' between runs and mislabels every map. Deriving them from the profile removes
#' the possibility.
#'
#' Clusters dominated by `oth` are ranked by `target_share`, the share of landed
#' **weight** in the three tunas, because that is what "non-target" means. Ranking
#' them by `oth_fraction` instead -- composition by *number* -- misnames a deep-set
#' tuna fishery that lands many small non-tuna as "non-target": in `sp_obshbf_lat`
#' the two `oth`-dominant clusters carry target shares of 0.132 and 0.676, and only
#' the first is not targeting tuna. A cluster above `target_max` is therefore called
#' "mixed", not "non-target".
#'
#' @param prof a [clusterProfile()] table.
#' @param target_max an `oth`-dominant cluster whose target share reaches this is
#'   called "mixed" rather than "non-target".
#' Best used once, at source: `cl <- factor(cl.K, labels = clusterLabels(prof))`.
#' Every downstream grid and plot then carries the fishery names automatically
#' and no `labels`/`prof` argument is needed anywhere.
#'
#' @return character vector, one label per cluster, in cluster order.
#' @family clustering utilities
#' @export
clusterLabels <- function(prof, target_max = 0.5) {
	nm  <- c(yft = "yellowfin", bet = "bigeye", alb = "albacore", skj = "skipjack")
	dom <- as.character(prof$dominant)
	out <- unname(nm[dom])
	isoth <- is.na(out)
	if (any(isoth)) {
		i <- which(isoth)
		key <- if (!is.null(prof$target_share)) prof$target_share[i]
			   else -prof$oth_fraction[i]           # fallback: number, not weight
		o  <- i[order(key)]                          # least tuna-directed first
		ts <- sort(key)
		lab <- character(length(o))
		nt <- ts < target_max
		lab[nt] <- if (sum(nt) == 1L) "non-target"
				   else c("pure non-target", "non-target 2",
						  paste("non-target", seq_len(max(0, sum(nt) - 2)) + 2))[seq_len(sum(nt))]
		lab[!nt] <- if (sum(!nt) == 1L) "mixed"
					else paste("mixed", seq_len(sum(!nt)))
		out[o] <- lab
	}

	## Labels MUST be unique. Two clusters can share a dominant species -- at
	## K = 6, sp_obshbf_lat returns two alb-dominant fisheries, alb 0.764 at
	## HBF 29 and alb 0.431 at HBF 20 -- and factor(cl, labels = ) MERGES levels
	## that share a label, silently turning six clusters into five. That is not
	## a cosmetic fault: the two fisheries are pooled in every downstream grid,
	## map and cross-tabulation, with no error raised anywhere.
	## Duplicates are suffixed in order of the dominant fraction, purest first.
	for (d in unique(out[duplicated(out)])) {
		i <- which(out == d)
		frac <- vapply(i, function(r) {
			cn <- paste0(dom[r], "_fraction")
			if (cn %in% names(prof)) prof[[cn]][r] else NA_real_
		}, 0)
		i <- i[order(frac, decreasing = TRUE)]
		out[i] <- paste(d, seq_along(i))
	}
	stopifnot(!anyDuplicated(out))
	out
}

plot.clusters <- function(cl.K, num_clusters = 3, reso=5){


	tag <- paste0(scenario,"_",num_clusters)
	figs.dir <- "figs/maps/"

	g <- clusterGrid(dat_clean, cl.K, species = "yft", res = 5)
	png(paste0(figs.dir,"cluster-dominance_",tag,".png"),1000,600)
	plotClusterDominance(g)
	dev.off()

	tag <- paste0(tag,"_",reso,"deg")
	g <- clusterGrid(dat_clean, cl.K, species = "yft", res = reso)
	png(paste0(figs.dir,"effort_",tag,".png"),1200,600)
	plotClusterMaps(g, what = "effort")
	dev.off()

	png(paste0(figs.dir,"yft-cpue_",tag,".png"),1200,600)
	plotClusterMaps(g, what = "cpue")
	dev.off()

	png(paste0(figs.dir,"yft-catch_",tag,".png"),1200,600)
	plotClusterMaps(g, what = "catch")
	dev.off()

	png(paste0(figs.dir,"hbf_",tag,".png"),1200,600)
	plotClusterMaps(g, what = "hbf")
	dev.off()	
}

flags.clusters <- function(cl.K, num_clusters = 3){

	tag <- paste0(scenario,"_",num_clusters)

	tb <- xtabs(E ~ flag + cl.K, dat_clean)
	
	f.out <- paste0(outdir,"/flags-",tag,".txt")
	write("#each fishery's fleet composition",f.out)
	write.table(round(prop.table(tb, 2), 3),f.out,quote=FALSE,
				col.names=T,row.names=T,sep="\t",append=T)		
	write("\n#each fleet's split across fisheries",f.out,append=T)
	write.table(round(prop.table(tb, 1), 3),f.out,quote=FALSE,
				col.names=T,row.names=T,sep="\t",append=T)		         
}



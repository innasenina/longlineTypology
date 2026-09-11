library(ggplot2)
for (f in list.files(pattern=".R", full.names = TRUE)) source(f)

scenarios <- c("yba_lat","sp","sp_lat","sp_hbf_lat","sp_obshbf_lat") 
Ks <- list(yba_lat=3,sp=c(4,5),sp_lat=5,sp_hbf_lat=5,sp_obshbf_lat=6)# <== K=5 is forced for sp_hbf_lat, otherwise with three components it returns K=4!
freq.tables <- list(yba_lat=table(rep(3,15)),sp=table(c(rep(4,8),rep(5,7))),
					sp_lat=table(rep(4,15)),sp_hbf_lat=table(rep(5,15)),
					sp_obshbf_lat=table(rep(5,15))) #<== K=5 forced for sp_hbf_lat and sp_obshbf_lat

#Run with 20,000 points:
#1) yba_lat
# pca n=2: var87, K=3
# pca n=3: var98, K=3

#2) sp
# pca n=2: var75, K=3
# pca n=3: var100, K=4

#3) sp_lat
# pca n=2: var69, K=3
# pca n=3: var89, K=4
# pca n=4: var100, K=5

#4) sp_obshbf_lat: 
# pca n=3: var78%, K=4
# pca n=4: var92%, K=6
# pca n=5: var100%, unstable, undefined - 7(2), max=10(13), ran with 100,000

f2use <- "old"				   

Ns <- 1

#What this script will execute. To run all scenarios:
#for (Ns in 1:5) source("../tests/test-gap.R")
scenario <- scenarios[Ns]
run.Gap  <- FALSE
run.Clus <- FALSE
run.diag <- TRUE

#options: define PCA by variance or number of components
set.pca.ncomp <- TRUE 
set.pca.var <- !set.pca.ncomp
write.RDS <- TRUE

# PCA and kmeans parameters
pca_variance_threshold <- 0.9   # Retain PCs explaining 70% of variance
sample_no              <- 20000 # Sample size per replicate for gap statistic
kmeans_set             <- 15    # Number of kmeans replicates
max_k                  <- 10    # Maximum number of clusters to test
iter_max               <- 1e6   # Maximum iterations for kmeans
nstart                 <- 30    # Number of random starts for kmeans
B_var                  <- 25    # Bootstrap replicates for gap statistic
nb_cores			   <- 5

# Cluster merging parameters
catch_threshold <- 3            # Merge clusters with < 3% of total catch

#Read pre-processed DATASET, it is already a data.frame
files <- list(old="~/WORK/Team/Romain/Fisheries-data/2026/LL_1x1_imp_raised_customLat.RDS",
			  new="~/WORK/Team/Romain/Fisheries-data/2026/LL_1x1_imp_lclim_raised_customLat.RDS")

dat_clean <- readRDS(files[[f2use]])

dat_clean <- dat_clean[which(dat_clean$date>=as.Date("1990-01-01")),]
if (length(grep("obshbf",scenario))>0)
	dat_clean <- dat_clean[which(dat_clean$imp_hbf==FALSE),]

if (f2use == "new") {
	colnames(dat_clean)[grep("ymd",colnames(dat_clean))] <- "date"
	#cut off 2023-2024 for consistency between old and new datasets
	dat_clean <- dat_clean[which(dat_clean$date<as.Date("2023-01-01")),]
}

## prepareEC() is idempotent, so calling it here costs nothing on data that
## has already been through it, and stops the silent failure where a missing
## shr simply drops target_share from the profile.

dat <- prepareEC(dat_clean)
shr <- targetShareW(dat)
if (is.null(shr))
	warning("no *_w columns: target_share will be absent from the profile")

dat <- setFractions(dat, scenario)

stat_list <- list()
stat_list$nrow_total <- nrow(dat)

# Step 1: Normalization and PCA
# Select and prepare variables
pca_select <- selectScenario(scenario)

cat("Variables selected for clustering:\n")
print(pca_select)

select_dat <- dat[, pca_select]

pca_ncomp  <- length(grep("_fraction$", pca_select, value = TRUE)) + ifelse(length(grep("_",scenario))>0,1,0) - 1  
if (set.pca.var)
	pca_res <- customPCA(pcaScale(select_dat, method = "zscore"), 
						 pca_variance_threshold, print_it = FALSE)
if (set.pca.ncomp)
	pca_res <- customPCA(pcaScale(select_dat, method = "zscore"), 
						 n_comp = pca_ncomp, print_it = FALSE)

#print retained variance by variable:
print(varRetained(pca_res, n_comp = pca_ncomp))

#get the variance from pca_res:
ev  <- pca_res$sdev^2
cum <- cumsum(ev) / sum(ev)
pca_var <- round(cum[pca_res$no_var]*100)	# % variance explained by the retained components

cat("PCA retained", pca_res$no_var, "components explaining ",
	pca_var, "% of variance\n\n")

pca_full <- pca_res$x[, 1:pca_res$no_var]


indices_90  <- which(as.integer(format(dat$date, "%Y")) >= 1990 & as.integer(format(dat$date, "%Y")) <= 2022)
#indices_90  <- which(as.integer(format(dat$date, "%Y")) >= 2023 & as.integer(format(dat$date, "%Y")) <= 2024)
pca_full_90 <- pca_full[indices_90, ]

# Step 2: Gap-statistic
if (run.Gap){
	set.seed(2)
	sample_indices <- lapply(1:kmeans_set, function(x) {
		sample(1:nrow(pca_full_90), size = sample_no, replace = TRUE)
	})

	cat("\nRunning", kmeans_set, "kmeans replicates on samples of",
		sample_no, "observations...\n")
		
	kmeans_res <- list()
	gap_stat   <- list()

	for (i in 1:kmeans_set) {
		message("i=",i)
		pca_sub <- pca_full_90[sample_indices[[i]], 1:pca_res$no_var]

		res <- customKmeans(pca_sub,
							max_k      = max_k,
							random_set = B_var,
							iter_max   = iter_max,
							nstart     = nstart,
							d.power    = 2,
							print_it   = FALSE)
							

		kmeans_res[[i]] <- res$kmeans
		gap_stat[[i]]   <- res$gap_stat
		plotGapSet(gap_stat, file = paste0("figs/gaps/gap_",scenario,
										   "_pca_n",pca_ncomp,"_var",pca_var,"_",sample_no,".png"))
	}

	# Determine most common K across replicates
	cluster_counts <- data.frame(
		cluster_n = sapply(kmeans_res, function(x) length(unique(x$cluster))),
		replicate = 1:kmeans_set
	)

	freq_table   <- table(cluster_counts$cluster_n)
	main_clust_n <- as.numeric(names(freq_table)[which.max(freq_table)])

	# Store as data frame for use in diagnostic report
	freq_df <- data.frame(
		n_clusters = as.numeric(names(freq_table)),
		freq       = as.numeric(freq_table)
	)

	cat("\nCluster number distribution across replicates:\n")
	print(freq_table)
	cat("\nSelected K =", main_clust_n, "clusters\n")
	stat_list$selected_k <- main_clust_n

	Ks.scenario <- as.integer(names(freq_table))
} else {

	Ks.scenario <- Ks[[scenario]] 
	freq_table  <- freq.tables[[scenario]]
}

out <- list(scenario = scenario, pca = pca_res, freq = freq_table, K = Ks.scenario,
			cluster = NULL, profile = NULL, row = NULL)
outdir <- "./outputs"
if (!dir.exists(outdir))
	dir.create(outdir, recursive = TRUE)

# Step 3: clustering
if (run.Clus){
	cat("Using", nrow(pca_full_90), "observations from 1990 onwards for clustering, 
		then assigning earlier records by the distance to cluster centroids\n")
	stat_list$nrow_clustering <- nrow(pca_full_90)

	for (K in Ks.scenario){

		#km.K <- kmeans(pca_full_90, centers = K, iter.max = iter_max, nstart = nstart)
		km.K <- kmeansProgress(pca_full_90, centers = K, iter.max = iter_max, nstart = nstart)

		cl.K <- assignClusters(as.data.frame(pca_full), as.data.frame(km.K$centers))


		z <- zeroTuna(dat)
		prof <- clusterProfile(dat, cl.K, share = shr, zero = z)

		#Label clusters for further plotting and diagnostics
		cl.K <- factor(cl.K, labels = clusterLabels(prof))

		#To write the cluster parameters into a summary file:
		out.row  <- scenarioRow(scenario, pca_res, freq_table, n_sub = sample_no)
		out.row$denominator <- attr(dat, "denominator")
		out.row$K_nontarget <- paste(nonTarget(prof), collapse = ",")
		out.row$K_dist <- if (run.Gap) out.row$K_dist else NA
		
		#To write the results to the out as RDS and write it optionally:
		out$cluster <- cl.K
		out$profile <- prof
		out$row<- out.row
		#write cluster's profiles into csv file
		write.csv(prof, file.path(outdir, paste0(scenario,"_",K,"_profile.csv")),
				  row.names = FALSE)
		f <- file.path(outdir, "summary.csv")
		write.table(out.row, f, sep = ",", row.names = FALSE,
					col.names = !file.exists(f), append = file.exists(f))
	
		if (write.RDS)
			saveRDS(out, file.path(outdir, paste0(scenario,"_",K,".RDS")))

		plot.clusters(cl.K,K,reso=1)
		flags.clusters(cl.K,K)
	}
}

if (!run.Clus & run.diag){
	for (K in Ks.scenario){
		fname <- paste0(outdir,"/",paste0(scenario,"_",K,".RDS"))
		cl  <- readRDS(fname)
		cl.K <- cl$cluster
		prof  <- cl$profile
		cl.lab <- factor(cl.K, labels = clusterLabels(prof))

		plot.clusters(cl.K,K,reso=1)
		flags.clusters(cl.K,K)
	}
}



# To compare clusters from different scenarios:
#adjRand(cl5_sp, cl5_sp_hbf_lat)
#nestedness(cl4_sp, cl5_sp)
  

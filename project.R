library(metacell)
library(stats)
library(bigmemory)
library(biganalytics)
library(ggplot2)
library(Seurat)
library(philentropy)
library(splatter)
library(DropletUtils)


## Functions

# Repeatedly runs k-means clustering on a KNN matrix and trims outliers each time 
# until it converges.
km_converge <- function(knn.mat, centers, max.iter=10) {
	filter.vec <- FALSE
	i <- 0
	while(any(!filter.vec) & (i < max.iter)) {
		size.co <- nrow(knn.mat) * .05 / centers
		knn.big <- as.big.matrix(knn.mat)
		km.clust <- bigkmeans(knn.big, centers=centers, iter.max=50, nstart=10, dist="euclid")
		cells <- rownames(knn.mat)
		assignments <- km.clust[1][[1]]
		i <- i + 1
		print(sprintf("%s iterations completed.", as.character(i)))
		# If any clusters remain with less than 25 entries, this process will repeat.
		filter.vec <- table(assignments) > size.co
		cell.filt <- filter.vec[assignments]
		knn.filt <- knn.mat[cell.filt, cell.filt]
		knn.mat <- knn.filt
	}
	return(list(cells, assignments))
}

# Boilerplate metacell preprocessing and KNN graph extraction.
preprocess_metacell <- function(assay, mat) {
	# Identify problematic genes
	nms <- c(rownames(mat@mat), rownames(mat@ignore_gmat))
	ig_genes <- c(grep("^IGJ", nms, v=T), grep("^IGH",nms,v=T), grep("^IGK", nms, v=T), 
				  grep("^IGL", nms, v=T))
	bad_genes <- unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),
		                  "NEAT1","TMSB4X", "TMSB10", ig_genes))
	# Boilerplaye metacell filtration KNN construction
	mcell_mat_ignore_genes(new_mat_id=assay, mat_id=assay, bad_genes, reverse=F) 
	mcell_mat_ignore_small_cells(assay, assay, 800)
	mcell_add_gene_stat(gstat_id=assay, mat_id=assay, force=T)
	mcell_gset_filter_varmean(gset_id="features", gstat_id=assay, T_vm=0.08, force_new=T)
	mcell_gset_filter_cov(gset_id="features", gstat_id=assay, T_tot=100, T_top3=2)
	mcell_add_cgraph_from_mat_bknn(mat_id=assay, gset_id="features", graph_id="graph",
								   K=50, dsamp=T)
	# Extract weighted balanced KNN Graph.
	knn.graph <- scdb_cgraph("graph")
	knn.from <- knn.graph@edges[[1]]
	knn.to <- knn.graph@edges[[2]]
	knn.weight <- knn.graph@edges[[3]]
	# This takes a long time.
	knn.mat <- unclass(xtabs(knn.weight ~ knn.from + knn.to))
	return(list(mat@mat, knn.mat))
}

# Constructs metacell expression vectors using the original MetaCell geometric mean formulation
make_metacells <- function(mat, partition) {
	metacells <- do.call("cbind", lapply(seq(max(partition)), function(part) {
		p.mat <- mat[, partition==part]
		if (class(p.mat)[1]=="numeric" | class(p.mat)[1]=="integer") {
			# Single-cell case.
			p.mat
		} else {
			# Cell number adjustment
			n.adj <- 1 / ncol(p.mat)
			# Geometric mean of each gene expression.
			apply(p.mat, 1, function(exp) {
				# Log-adjusted
				logxp <- sum(sapply(exp, function(c.exp) {log(1 + c.exp)}))
				numerator <- exp(n.adj * logxp) - 1
				# Adjusts for average size factor of composite cells
				denominator <- n.adj * sum(p.mat)
				numerator / denominator
			})
		}
	}))
	return(metacells)
}

# Counts non-zero entries of each column of a matrix, returns as vector.
mat_sparsity <- function(mat) {
	apply(mat, 2, function(col){sum(col > 0)})
}

# Balancing metacells per cluster
make_balanced_metacells <- function(seurat.object, knn, count.mat) {
	clusters <- seurat.object$seurat_clusters
	split.data <- lapply(levels(clusters), function(cluster) {
		# SCDB preprocessing
		cell.mask <- (clusters==cluster)
		split.knn <- knn[cell.mask, cell.mask]
		# centers = 25
		km.data.10 <- km_converge(split.knn, centers=10, max.iter=10)
		km.cells.10 <- km.data.10[[1]]
		km.ass.10 <- km.data.10[[2]]
		mat.filt.10 <- count.mat[, km.cells.10]
		metacells <- make_metacells(mat.filt.10, km.ass.10)
		colnames(metacells) <- paste0("c", cluster, "m", seq(10))
		list(metacells, km.ass.10)
	})
	split.mc <- lapply(split.data, function(dat) {dat[[1]]})
	split.ass <- lapply(split.data, function(dat) {dat[[2]]})
	return(list(split.mc, split.ass))
}

# Euclidean distance formula
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# Creates a density vector of distances between points in two given matrices.
mat_dist <- function(mat.1, mat.2) {
	dist <- unlist(apply(mat.1, 1, function(row.1){
	unlist(apply(mat.2, 1, function(row.2) {euclidean(row.1, row.2)}))
	}))
	freq <- hist(dist, breaks=seq(min(dist), max(dist), length.out=11))$counts
	freq / sum(freq)
}

# Converts between cluster labels in two datasets
cluster_convert <- function(clust.8k, conv.6k=c(1, 3, 5, 4, 2, 7, 6, 8, 9)) {
	clust.6k <- conv.6k[clust.8k + 1] - 1
	return(clust.6k)
}

# Extracts mean gene expression from Splatter simulation
gene_mean <- function(sim.g.meta, gene, batch, ct) {
	base.mean <- sim.g.meta[gene, "GeneMean"]
	batch.id <- paste0("BatchFacBatch", batch)
	batch.eff <- sim.g.meta[gene, batch.id]
	ct.id <- paste0("DEFacGroup", ct)
	ct.eff <- sim.g.meta[gene, ct.id]
	adjusted.mean <- base.mean * batch.eff * ct.eff / 50000
	return(adjusted.mean)
}


# Constructs metacells with a given number of centers either from an entire dataset or from 
# a predefined partition of that dataset.
make_metacells_generic <- function(count.mat, centers, partition=NULL) {
	# Write count.mat as a 10X object to load into SCDB
	sparse.count.mat <- as(count.mat, "dgCMatrix")
	write10xCounts("temp/", sparse.count.mat, overwrite=TRUE)
	assay <- "t"
	scdb_init("temp/", force_reinit=T)
	mcell_import_scmat_10x(assay, base_dir="temp/")
	mat <- scdb_mat(assay)
	pre.data <- preprocess_metacell(assay, mat)
	pre.mat <- pre.data[[1]]
	knn <- pre.data[[2]]
	if (is.null(partition)) {
		km.data <- km_converge(knn, centers=centers)
		km.cells <- km.data[[1]]
		km.ass <- km.data[[2]]
		mat.filt <- pre.mat[, km.cells]
		metacells <- make_metacells(mat.filt, km.ass)
		return(list(metacells, km.ass))
	} else {
		split.data <- lapply(levels(partition), function(part) {
			# SCDB preprocessing
			cell.mask <- (partition==part)
			split.knn <- knn[cell.mask, cell.mask]
			# centers = 25
			km.data <- km_converge(split.knn, centers=centers, max.iter=10)
			km.cells <- km.data[[1]]
			km.ass <- km.data[[2]]
			mat.filt <- count.mat[, km.cells]
			metacells <- make_metacells(mat.filt, km.ass)
			colnames(metacells) <- paste0(part, "m", seq(10))
			list(metacells, km.ass)
		})
		split.mc <- lapply(split.data, function(dat) {dat[[1]]})
		split.ass <- lapply(split.data, function(dat) {dat[[2]]})
		return(list(split.mc, split.ass))
	}
}

# Gets ground truth metacell expression levels by weighted mean of expression levels from Splatter.
get_metacell_gts <- function(mean.frame, cell.types, assignments) {
	do.call("cbind", lapply(levels(as.factor(assignments)), function(ass) {
		cell.mask <- (assignments==ass)
		ct.count <- as.vector(table(cell.types[cell.mask]))
		ct.prop <- ct.count / sum(ct.count)
		a.means <- sweep(mean.frame, MARGIN=2, ct.prop, `*`)
		gt.vec <- rowSums(a.means)
	}))
}

# Calculates average deviation from ground truth of each metacell.
get_metacell_res <- function(mean.frame, cell.types, mc.data) {
	mc <- mc.data[[1]]
	assignments <- mc.data[[2]]
	mc.gt <- get_metacell_gts(mean.frame, cell.types, assignments)
	mc.dev <- mc - mc.gt
	mc.res <- abs(colMeans(mc.dev) + (2e-5))
}

# Set file paths
source.dir <- "/scratch/tmp/sereno/raphael/cos597final/sourcefiles"
fig.dir <- "/scratch/tmp/sereno/raphael/cos597final/figs"
source.6k <- file.path(source.dir, "pbmc_6")
source.8k <- file.path(source.dir, "pbmc_8")
data.dir <- "/scratch/tmp/sereno/raphael/cos597final/data"
data.6k <- file.path(data.dir, "pbmc_6")
data.8k <- file.path(data.dir, "pbmc_8")

# Process 6k
# assay <- "6k"
scdb_init(data.6k, force_reinit=T)
mcell_import_scmat_10x("6k", base_dir=source.6k)
mat <- scdb_mat(assay)
pre.data <- preprocess_metacell(assay, mat)
pre.mat <- pre.data[[1]]
knn <- pre.data[[2]]
# centers = 25
km.data.25 <- km_converge(knn, centers=25)
km.cells.25 <- km.data.25[[1]]
km.ass.25 <- km.data.25[[2]]
mat.filt.25 <- pre.mat[, km.cells.25]
metacells.25 <- make_metacells(mat.filt.25, km.ass.25)
# centers = 50
km.data.50 <- km_converge(knn, centers=50)
km.cells.50 <- km.data.50[[1]]
km.ass.50 <- km.data.50[[2]]
mat.filt.50 <- pre.mat[, km.cells.50]
metacells.50 <- make_metacells(mat.filt.50, km.ass.50)
# centers = 100
km.data.100 <- km_converge(knn, centers=100)
km.cells.100 <- km.data.100[[1]]
km.ass.100 <- km.data.100[[2]]
mat.filt.100 <- pre.mat[, km.cells.100]
metacells.100 <- make_metacells(mat.filt.100, km.ass.100)


scdb_init(data.6k, force_reinit=T)
mcell_import_scmat_10x("6k", base_dir=source.6k)
mat.6k <- scdb_mat("6k")
pre.data.6k <- preprocess_metacell("6k", mat.6k)
pre.mat.6k <- pre.data.6k[[1]]
knn.6k <- pre.data.6k[[2]]

scdb_init(data.8k, force_reinit=T)
mcell_import_scmat_10x("8k", base_dir=source.8k)
mat.8k <- scdb_mat("8k")
pre.data.8k <- preprocess_metacell("8k", mat.8k)
pre.mat.8k <- pre.data.8k[[1]]
knn.8k <- pre.data.8k[[2]]


# Seurat PCA and Louvrain on single cells
seur.data.6k <- Read10X(source.6k)
seur.data.6k.filt <- seur.data.6k[, match(colnames(pre.mat.6k), colnames(seur.data.6k))]
seur.6k <- CreateSeuratObject(seur.data.6k.filt)
seur.6k <- NormalizeData(seur.6k)
seur.6k <- FindVariableFeatures(seur.6k, nfeatures=2000)
seur.6k <- ScaleData(seur.6k)
seur.6k <- RunPCA(seur.6k, features=VariableFeatures(object=seur.6k))
seur.6k <- FindNeighbors(seur.6k, dims=1:15)
seur.6k <- FindClusters(seur.6k, resolution=0.4)
s.m.6k.dat <- make_balanced_metacells(seur.6k, knn.6k, pre.mat.6k)
s.m.6k <- s.m.6k.dat[[1]]
s.ass.6k <- s.m.6k.dat[[2]]
s.agg.6k <- unlist(lapply(s.ass.6k, function(ass) {
	as.vector(table(ass))
}))
# Post-metacell processing
pca.6k <- Embeddings(seur.6k)[, 1:20]
mc.data.6k <- do.call("cbind", s.m.6k)
mc.seur.6k <- CreateSeuratObject(mc.data.6k)
mc.seur.6k <- FindVariableFeatures(mc.seur.6k, nfeatures=2000)
mc.seur.6k <- ScaleData(mc.seur.6k)
mc.seur.6k <- RunPCA(mc.seur.6k, features=VariableFeatures(object=mc.seur.6k))
mc.pca.6k <- Embeddings(mc.seur.6k)[, 1:20]


seur.data.8k <- Read10X(source.8k)
# seur.data.8k.filt <- seur.data.8k[, match(colnames(mat.8k@mat), colnames(seur.data.8k))]
seur.8k <- CreateSeuratObject(seur.data.8k)
seur.8k <- NormalizeData(seur.8k)
seur.8k <- FindVariableFeatures(seur.8k, nfeatures=2000)
seur.8k <- ScaleData(seur.8k)
seur.8k <- RunPCA(seur.8k, features=VariableFeatures(object=seur.8k))
seur.8k <- FindNeighbors(seur.8k, dims=1:10)
seur.8k <- FindClusters(seur.8k, resolution=0.2)
s.m.8k.dat <- make_balanced_metacells(seur.8k, knn.8k, pre.mat.8k)
s.m.8k <- s.m.8k.dat[[1]]
s.ass.8k <- s.m.8k.dat[[2]]
s.agg.8k <- unlist(lapply(s.ass.8k, function(ass) {
	as.vector(table(ass))
}))
# Post-metacell processing.
pca.8k <- Embeddings(seur.8k)[, 1:20]
mc.data.8k <- do.call("cbind", s.m.8k)
mc.seur.8k <- CreateSeuratObject(mc.data.8k)
mc.seur.8k <- FindVariableFeatures(mc.seur.8k, nfeatures=2000)
mc.seur.8k <- ScaleData(mc.seur.8k)
mc.seur.8k <- RunPCA(mc.seur.8k, features=VariableFeatures(object=mc.seur.8k))
mc.pca.8k <- Embeddings(mc.seur.8k)[, 1:20]
# unique(seur.8k$seurat_clusters)


## Splatter Simulation

params <- newSplatParams()
cluster.probs <- c(0.22, 0.19, 0.17, 0.14, 0.12, 0.08, 0.05, 0.02, 0.01)
sim <- splatSimulate(batchCells=c(5000, 5000), group.prob=cluster.probs, method="groups")
sim.counts <- counts(sim)
sim.c.meta <- colData(sim)
sim.ncounts <- t(counts(sim)) / sim.c.meta$ExpLibSize
sim.g.meta <- rowData(sim)
# List of gene x ct frames with each list index for a batch. Lookup table for expected mean.
mean.frames <- lapply(seq(2), function(b) {
	do.call("rbind", lapply(seq(nrow(sim.g.meta)), function(g) {
		b.frame <- do.call("cbind", lapply(seq(9), function(c) {
				gene_mean(sim.g.meta, g, b, c)
		}))
		colnames(b.frame) <- paste0("B", b, "C", seq(9))
		b.frame
	}))
})
# Calculate and plot residuals
cell.residuals <- sapply(seq(nrow(sim.c.meta)), function(cell) {
	cell.meta <- sim.c.meta[cell, ]
	cell.batch <- ifelse(cell.meta$Batch=="Batch1", 1, 2)
	cell.ct <- as.numeric(gsub("Group", "", cell.meta$Group))
	mean.vec <- mean.frames[[cell.batch]][, cell.ct]
	true.vec <- sim.ncounts[cell, ]
	residual <- mean(true.vec - mean.vec)
})
res.frame <- data.frame(residual=cell.residuals, Batch=rep(c("Batch 1", "Batch 2"), each=5000))
sci.res <- ggplot(data=res.frame, aes(x=cell.residuals)) + 
    geom_density() +
    facet_wrap(~Batch) + 
    # scale_fill_brewer(palette="Dark2") + 
    ggtitle("Deviation from Ground Truth Expression by Batch, Single-Cells") + 
  	xlab("Deviation from Ground Truth Mean Expression") + 
    ylab("Frequency") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    	panel.border=element_rect(colour="grey89", fill=NA, size=0.5))
sci.res.path <- file.path(fig.dir, "sci_res.png")
ggsave(filename=sci.res.path, plot=sci.res, width=10, height=5, units="in", dpi=400)
# Split data by batch
b1.mask <- sim.c.meta$Batch=="Batch1"
b2.mask <- sim.c.meta$Batch=="Batch2"
counts.1 <- sim.counts[, b1.mask]
counts.2 <- sim.counts[, b2.mask]
ncounts.1 <- sim.ncounts[b1.mask, ]
ncounts.2 <- sim.ncounts[b2.mask, ]
cm.1 <- sim.c.meta[b1.mask, ]
cm.2 <- sim.c.meta[b2.mask, ]

split.data.1 <- make_metacells_generic(counts.1, centers=10, partition=cm.1$Group)
mc.res.1 <- unlist(lapply(seq(9), function(g) {
	mc <- split.data.1[[1]][[g]]
	gt <- mean.frames[[1]][, g]
	# Adjustment for baseline shift of Splatter's algorithm
	apply(mc, 2, function(exp) {abs(mean(exp - gt) + (2e-5))})
}))
# Outlier clipping
mc.res.1[mc.res.1 > 7e-6] <- 7e-6
mc.agg.1 <- unlist(lapply(split.data.1[[2]], function(ass) {as.vector(table(ass))}))
split.data.2 <- make_metacells_generic(counts.2, centers=10, partition=cm.2$Group)
mc.res.2 <- unlist(lapply(seq(9), function(g) {
	mc <- split.data.2[[1]][[g]]
	gt <- mean.frames[[2]][, g]
	apply(mc, 2, function(exp) {abs(mean(exp - gt) + (2e-5))})
}))
mc.res.2[mc.res.2 > 7e-6] <- 7e-6
mc.agg.2 <- unlist(lapply(split.data.2[[2]], function(ass) {as.vector(table(ass))}))
res.frame <- data.frame(residual=c(mc.res.1, mc.res.2), size=c(mc.agg.1, mc.agg.2), 
	Batch=rep(c("Batch 1", "Batch 2"), each=90), 
	Cluster=rep(rep(paste0("Cell Type ", seq(9)), each=10), 2))
sim.spar <- ggplot(data=res.frame, aes(x=size, y=residual)) + 
    geom_point(size=2, aes(color=Cluster)) +
    facet_wrap(~Batch) + 
    scale_fill_brewer(palette="Dark2") + 
    ggtitle("Deviation From Ground Truth Expression by Metacell Size") + 
  	xlab("Number of Cells in Aggregation") + 
    ylab("Average Deviation from Ground Truth Expression") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    	panel.border=element_rect(colour="grey89", fill=NA, size=0.5))
sim.spar.path <- file.path(fig.dir, "sim_residual_10.png")
ggsave(filename=sim.spar.path, plot=sim.spar, width=11, height=5, units="in", dpi=400)


mc.data.25.1 <- make_metacells_generic(counts.1, centers=25)
res.25.1 <- get_metacell_res(mean.frames[[1]], cm.1$Group, mc.data.25.1)
agg.25.1 <- as.vector(table(mc.data.25.1[[2]]))
mc.data.25.2 <- make_metacells_generic(counts.2, centers=25)
res.25.2 <- get_metacell_res(mean.frames[[2]], cm.2$Group, mc.data.25.2)
agg.25.2 <- as.vector(table(mc.data.25.2[[2]]))
# 50 Centers
mc.data.50.1 <- make_metacells_generic(counts.1, centers=50)
res.50.1 <- get_metacell_res(mean.frames[[1]], cm.1$Group, mc.data.50.1)
res.50.1[res.50.1 > 7e-6] <- 7e-6
agg.50.1 <- as.vector(table(mc.data.50.1[[2]]))
mc.data.50.2 <- make_metacells_generic(counts.2, centers=50)
res.50.2 <- get_metacell_res(mean.frames[[2]], cm.2$Group, mc.data.50.2)
res.50.2[res.50.2 > 7e-6] <- 7e-6
agg.50.2 <- as.vector(table(mc.data.50.2[[2]]))

size.vec <- c(agg.25.1, agg.50.1, mc.agg.1, agg.25.2, agg.50.2, mc.agg.2)
res.vec <- c(res.25.1, res.50.1, mc.res.1, res.25.2, res.50.2, mc.res.2)
id.char <- c(rep("25 Metacells", 25), rep("50 Metacells", 50), 
	rep("10 Metacells Per\nCell Type (Total 90)", 90))
id.vec <- factor(rep(id.char, 2),
    levels=c("25 Metacells", "50 Metacells", "10 Metacells Per\nCell Type (Total 90)"))
batch.vec <- rep(c("Batch 1", "Batch 2"), each=165)
res.k.frame <- data.frame(size=size.vec, residual=res.vec, Aggregation=id.vec, Batch=batch.vec)
res.plot <- ggplot(data=res.k.frame, aes(x=size, y=residual)) + 
    geom_point(size=2, aes(color=Aggregation)) +
    facet_wrap(~Batch) + 
    scale_fill_brewer(palette="Dark2") + 
    ggtitle("Deviation From Ground Truth Expression by Metacell Size") + 
  	xlab("Number of Cells in Aggregation") + 
    ylab("Average Deviation from Ground Truth Expression") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    	panel.border=element_rect(colour="grey89", fill=NA, size=0.5))
res.path <- file.path(fig.dir, "batch_residual.png")
ggsave(filename=res.path, plot=res.plot, width=11, height=5, units="in", dpi=400)


## Didn't work.
# k-l divergence of distances between single cells vs metacells for two clusters
clust.6k <- seur.6k$seurat_clusters
clust.8k <- seur.8k$seurat_clusters
lapply(seq(0, 8), function(i) {
	do.call("rbind", lapply(seq(0, 8), function(j) {
		c.1 <- cluster_convert(i)
		c.2 <- cluster_convert(j)
		c.3 <- i
		c.4 <- j
		# Single-cell KL
		sc.1 <- seur.data.6k.filt[, clust.6k==c.1]
		sc.2 <- seur.data.6k.filt[, clust.6k==c.2]
		sc.dist.1 <- mat_dist(sc.1[, sample(ncol(sc.1), 10)], sc.2[, sample(ncol(sc.2), 10)])
		sc.3 <- seur.data.8k.filt[, clust.8k==c.1]
		sc.4 <- seur.data.8k.filt[, clust.8k==c.2]
		sc.dist.2 <- mat_dist(sc.3[, sample(ncol(sc.3), 10)], sc.4[, sample(ncol(sc.4), 10)])
		sc.kl <- KL(rbind(sc.dist.1, sc.dist.2), unit="log")
		# Metacell KL
		seq.1 <- seq((c.1*10 + 1), (c.1*10+10))
		seq.2 <- seq((c.2*10 + 1), (c.2*10+10))
		mc.1 <- mc.data.6k[, seq.1]
		mc.2 <- mc.data.6k[, seq.2]
		mc.dist.1 <- mat_dist(mc.1, mc.2)
		seq.3 <- seq((c.3*10 + 1), (c.3*10+10))
		seq.4 <- seq((c.4*10 + 1), (c.4*10+10))
		mc.3 <- mc.data.8k[, seq.1]
		mc.4 <- mc.data.8k[, seq.2]
		mc.dist.2 <- mat_dist(mc.3, mc.4)
		mc.kl <- KL(rbind(mc.dist.1, mc.dist.2), unit="log")
		data.frame(kl=c(sc.kl, mc.kl), id=c("sc", "mc"))
	}))
})
# Labeling
markers.6k <- FindAllMarkers(seur.6k)
unlist(sapply(seq(0, 8), function(n) {
	m <- markers.6k[markers.6k$cluster==n, ]
	m[1, ]
})["gene", ])
markers.8k <- FindAllMarkers(seur.8k)
unlist(sapply(seq(0, 8), function(n) {
	m <- markers.8k[markers.8k$cluster==n, ]
	m[1, ]
})["gene", ])
# 1, 2, 3, 4, 5, 6, 7, 8, 9
# "RPL32"  "IL32"   "S100A8" "CCL5"   "CD79A"  "FCGR3A" "GNLY"   "FCER1A" "PF4"
# 1, 3, 5, 4, 2, 7, 6, 8, 9
# "RPL32"  "S100A8" "CD79A"  "CCL5"   "IL32"   "GNLY"   "FCGR3A" "FCER1A" "PF4"


## PLOTS
# Sparsity of each set of cells or metacells
sparsity.vec <- c(mat_sparsity(pre.mat), mat_sparsity(metacells.25), mat_sparsity(metacells.50), 
	mat_sparsity(mc.data.8k))
# Number of cells aggregated into each one.
agg.vec <- c(rep(1, ncol(pre.mat)), as.vector(table(km.ass.25)), as.vector(table(km.ass.50)), 
	s.agg.8k)
# Original datasets of each.
dataset <- factor(c(rep("Single-Cell", ncol(pre.mat)), rep("25 Metacells", 25), 
	rep("50 Metacells", 50), 
	rep("10 Metacells Per\nCell Type (Total 90)", 90)), 
	levels=c("25 Metacells", "50 Metacells", "10 Metacells Per\nCell Type (Total 90)",
		"Single-Cell"))
spar.frame <- data.frame(sparsity=sparsity.vec, count=agg.vec, Aggregation=dataset)
spar.scatter <- ggplot(data=spar.frame, aes(x=count, y=sparsity)) + 
    geom_point(size=2, aes(color=Aggregation)) +
    scale_fill_brewer(palette="Dark2") + 
    ggtitle("Nonzero Entries by Cell Aggregation, PBMC 8K") + 
  	xlab("Number of Cells in Aggregation") + 
    ylab("Nonzero Entries") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    	panel.border=element_rect(colour="grey89", fill=NA, size=0.5))
spar.path <- file.path(fig.dir, "sparsity_8k.png")
ggsave(filename=spar.path, plot=spar.scatter, width=8, height=6, units="in", dpi=400)
# Plot just by-cluster metacells
ct.sparsity.vec <- mat_sparsity(mc.data.8k)
ct.agg.vec <- s.agg.8k
cell.types <- paste("Cell Type", seq(9))
cluster <- rep(cell.types, each=10)
ct.frame <- data.frame(sparsity=ct.sparsity.vec, count=ct.agg.vec, Cell_Type=cluster)
ct.scatter <- ggplot(data=ct.frame, aes(x=count, y=sparsity)) + 
    geom_point(size=2, aes(color=Cell_Type)) +
    scale_fill_brewer(palette="Dark2") + 
    ggtitle("Nonzero Entries by Cell Type, PBMC 8K") + 
  	xlab("Number of Cells in Aggregation") + 
    ylab("Nonzero Entries") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    	panel.border=element_rect(colour="grey89", fill=NA, size=0.5))
ct.path <- file.path(fig.dir, "ct_sparsity_8k.png")
ggsave(filename=ct.path, plot=ct.scatter, width=8, height=6, units="in", dpi=400)

library(ArchR)
set.seed(1)

# Configure
addArchRThreads(threads = 15) 
addArchRGenome('hg19')


# ################################################################################################
# Arrow files and project 

# Input files
data_dir = <Directory containing the ATAC fragments file>
setwd(sprintf("%s/ArchR", data_dir))
inputFiles <- c(sprintf("%s/atac_pbmc_10k_nextgem_fragments.tsv.gz", data_dir)
              )
names(inputFiles) <- c('pbmc_10k_atac'
                       )

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 1, #Dont set this too high because you can always increase later
  filterFrags =3000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
  removeFilteredCells = TRUE
)


# Create project
proj_name <- "pbmc_10x_atac"
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = proj_name,
  copyArrows = FALSE
)


# ################################################################################################
# PReprocesing

# SVD, Clustering, UMAP
res <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                       name = "IterativeLSI", scaleDims=FALSE, force=TRUE)#, varFeatures=100000)
proj <- res[[1]]
var_features <- res[[2]]

# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)



# Peaks 
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addGroupCoverages(proj, maxFragmentLength=147)
proj <- addReproduciblePeakSet(proj)
# Counts
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)

# Save 
proj <- saveArchRProject(ArchRProj = proj)



# ################################################################################################

# Export
dir.create(sprintf("%s/export", proj_name))
write.csv(getReducedDims(proj), sprintf('%s/export/svd.csv', proj_name), quote=FALSE)
write.csv(getCellColData(proj), sprintf('%s/export/cell_metadata.csv', proj_name), quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
scores <- as.matrix(scores)
rownames(scores) <- rowData(gene.scores)$name
write.csv(scores, sprintf('%s/export/gene_scores.csv', proj_name), quote=FALSE)



# Peak counts
peaks <- getPeakSet(proj)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')

# Reorder peaks 
# Chromosome order
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order)
    reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
reordered_features <- Reduce("c", reordered_features)    

# Export counts
dir.create(sprintf("%s/export/peak_counts", proj_name))
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, sprintf('%s/export/peak_counts/counts.mtx', proj_name))
write.csv(colnames(peak.counts), sprintf('%s/export/peak_counts/cells.csv', proj_name), quote=FALSE)
names(reordered_features) <- sprintf("Peak%d", 1:length(reordered_features))
write.csv(as.data.frame(reordered_features), sprintf('%s/export/peak_counts/peaks.csv', proj_name), quote=FALSE)



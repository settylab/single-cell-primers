library(ArchR)
library(parallel)
set.seed(1)
#system("ml MACS2/2.2.6-foss-2019b-Python-3.7.4")
# Configure
addArchRThreads(threads = 12)
addArchRGenome('hg38')


# ################################################################################################
# Arrow files and project 

# Input files
#data_dir = '/fh/fast/setty_m/user/cjordan2/repositories/single-cell-primers/data/ArchR'
#output_dir = '/fh/fast/setty_m/user/cjordan2/repositories/single-cell-primers/tcell-multiome-data'
#setwd(data_dir)

# Input files
data_dir = '/fh/fast/setty_m/grp/lab-datasets/bonemarrow-tcell-dep-multiome'
output_dir = '/fh/fast/setty_m/user/cjordan2/repositories/single-cell-primers'
#dir.create(sprintf("%s/ArchR", output_dir))
inputFiles <- c(
    sprintf("%s/cr-arc-results/rep1/atac_fragments.tsv.gz", data_dir),
    sprintf("%s/cr-arc-results/rep2/atac_fragments.tsv.gz", data_dir)
              )
sample <- c('rep1', 'rep2')
names(inputFiles) <- sample

# Subset of cells detcrmined in RNA
# Multiome
multiome_path =  sprintf("%s/data", output_dir)
valid_barcodes = list()
rep1 = read.csv(sprintf("%s/tcell-dep_atac_cells_rep1.csv", multiome_path), stringsAsFactors = FALSE)[,1]
valid_barcodes[[1]] = rep1
rep2 = read.csv(sprintf("%s/tcell-dep_atac_cells_rep2.csv", multiome_path), stringsAsFactors = FALSE)[,1]
valid_barcodes[[2]] = rep2
#valid_barcodes = getValidBarcodes(c(sprintf("%s/rep1_barcodes.csv", multiome_path), sprintf("%s/rep2_barcodes.csv", multiome_path)), sampleNames = c("rep1", "rep2"))
setwd(sprintf("%s/ArchR", output_dir))
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  #minTSS = 1, 
  #minFrags = 500, 
  validBarcodes = valid_barcodes,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
  force = FALSE
)
#proj = loadArchRProject(data_dir+ 'tcell_dep')
#proj_name <- "temp"
#proj <- ArchRProject(
#  ArrowFiles = ArrowFiles, 
#  outputDirectory = proj_name,
#  copyArrows = FALSE
#)


# Create Arrow files 
# Note that the TSS and Frags filter might result in some cells not being included. 
# Set these to 0 if you would like all cells to included.

#multiome_cells <- intersect(multiome_cells, getCellNames(proj))

#proj_name <- "tcell_atac"
#proj <- ArchRProject(
#  ArrowFiles = ArrowFiles, 
#  outputDirectory = proj_name,
#  copyArrows = FALSE,
#)


proj <- loadArchRProject(sprintf("%s/ArchR/tcell_atac", output_dir))
# ################################################################################################
# PReprocesing

# SVD, Clustering, UMAP\
fun <- function(){
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                       name = "IterativeLSI", scaleDims=FALSE, force=FALSE, varFeatures=100000)
var_features <- proj@reducedDims[["IterativeLSI"]]$LSIFeatures

# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)

proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', threads = 1, force = TRUE)



# Peaks 
proj <- addClusters(input = proj, maxFragmentLength = 147, reducedDims = "IterativeLSI", force = FALSE)
proj <- addGroupCoverages(proj, force = FALSE, threads = 1)
proj <- addReproduciblePeakSet(proj, force = FALSE)
# Counts
proj <- addPeakMatrix(proj,ceiling=10^9, maxFragmentLength = 147, force = TRUE)


#imputation
proj <- addImputeWeights(proj)

# ChromVAR

#system(sprintf("rm %s/tcell_dep/Background-Peaks.rds", data_dir))
proj <- addBgdPeaks(proj,method = 'chromVAR', force = TRUE)
proj <- addMotifAnnotations(ArchRProj= proj, motifSet="cisbp", name="Motif", force=TRUE)
proj <- addDeviationsMatrix(ArchRProj= proj, peakAnnotation="Motif", threads=1, force=TRUE)



## chromvar



# Motif ranking
var_deviations <- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE)
markerMotifs <- sprintf("z:%s", var_deviations$name)
warnings()
# Motif scores
motifScores <- ArchR:::.getMatrixValues(proj, name = markerMotifs, matrixName = 'MotifMatrix', log2Norm = FALSE)
# Save 
proj <- saveArchRProject(ArchRProj = proj)



}
# ################################################################################################
proj_name = sprintf("%s/ArchR/tcell_atac", output_dir)
# Export
dir.create(sprintf("%s/export", proj_name))
write.csv(motifScores, sprintf('%s/export/motif_scores.csv', proj_name), quote=FALSE)


message("successfully ChromVAR exported")
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


#imputed gene_scores
imputed.weights <-getImputeWeights(proj)
write.csv(imputed.weights, sprintf('%s/export/imputed_weights.csv', proj_name), quote=FALSE)
imputed_scores <- imputeMatrix(mat = assay(gene.scores), imputeWeights = imputed.weights)
write.csv(imputed_scores, sprintf('%s/export/imputed_scores.csv', proj_name), quote=FALSE)

scores <-assays(gene.scores)
# chromVar
var_dev<- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE)
write.csv(var_dev, sprintf("%s/export/deviations.csv", proj_name), quote = FALSE)
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



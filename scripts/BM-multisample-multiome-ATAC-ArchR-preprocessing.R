library(ArchR)
set.seed(1)

# Configure
addArchRThreads(threads = 1) 
addArchRGenome('hg38')


# ################################################################################################
# Arrow files and project 

# Input files
# data_dir = <Directory containing the ATAC fragments file>
data_dir = '/fh/fast/setty_m/user/msetty/repositories/single-cell-primers/notebooks/multi-sample-multiome/'
inputFiles <- c(
    sprintf("%s/rep1/atac_fragments.tsv.gz", data_dir),
    sprintf("%s/rep2/atac_fragments.tsv.gz", data_dir)
              )
samples <- c('bm_multiome_rep1', 'bm_multiome_rep2')
names(inputFiles) <- c(samples)

# Subset of cells determined in RNA
# Multiome
# multiome_path = <Directory where multiome results were exported>
multiome_path = 'multi-sample-multiome/'
multiome_cells = read.csv(sprintf("%s/bm_multiome_cells.csv", multiome_path), stringsAsFactors=FALSE)[,2]
valid_barcodes = list()
for (sample in samples){
    valid_barcodes[[sample]] = gsub('.*#', '', multiome_cells[grep(sample, multiome_cells)])
}

dir.create(sprintf("%s/ArchR", data_dir))
setwd(sprintf("%s/ArchR", data_dir))


# Create Arrow files 
# Note that the TSS and Frags filter might result in some cells not being included. 
# Set these to 0 if you would like all cells to included.
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 1, 
  minFrags = 500, 
  validBarcodes = valid_barcodes,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
)


# Create project
proj_name <- "bm_multiome_atac"
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = proj_name,
  copyArrows = FALSE
)


# ################################################################################################
# Preprocesing

# Gene scores
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix')

# IterativeLSI
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# Peaks 
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addGroupCoverages(proj, maxFragmentLength=147)
proj <- addReproduciblePeakSet(proj)
# Counts
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)

# UMAPs
proj <- addUMAP(proj)

# Save 
proj <- saveArchRProject(ArchRProj = proj)



# ################################################################################################

# Export
dir.create(sprintf("%s/export", proj_name))
write.csv(getReducedDims(proj), sprintf('%s/export/svd.csv', proj_name), quote=FALSE)
write.csv(getCellColData(proj), sprintf('%s/export/cell_metadata.csv', proj_name), quote=FALSE)
write.csv(getEmbedding(proj), sprintf('%s/export/umap.csv', proj_name), quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
dir.create(sprintf("%s/export/gene_scores", proj_name))
writeMM(scores, sprintf('%s/export/gene_scores/scores.mtx', proj_name))
write.csv(colnames(scores), sprintf('%s/export/gene_scores/cells.csv', proj_name), quote=FALSE)
write.csv(rowData(gene.scores)$name, sprintf('%s/export/gene_scores/genes.csv', proj_name), quote=FALSE)



# Peak counts
peaks <- getPeakSet(proj)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')

# Reorder peaks 
# Chromosome order [This mess is necessary since the peaks get sorted by lexicographical order]
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





# ################################################################################################

# chromVAR scores
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
# Motif scores
motifScores <-   getMatrixFromProject(proj, 'MotifMatrix')
scores <- as.data.frame(assays(motifScores)[['z']])
write.csv(scores, quote=FALSE,
         sprintf("%s/export/chromvar_motif_scores.csv", proj_name))

# Save the project
proj <- saveArchRProject(ArchRProj = proj)

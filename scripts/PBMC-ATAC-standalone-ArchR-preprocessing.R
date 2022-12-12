library(ArchR)
set.seed(1)

# Number of threads and gnome
num_threads <- 1
genome <- 'hg19'

# Configure
addArchRThreads(threads = num_threads) 
addArchRGenome(genome)


# ################################################################################################
# Arrow files and project 

# Input files
data_dir = "/fh/fast/setty_m" #redirect tocorrect dir
setwd(sprintf("%s/ArchR", data_dir))
inputFiles <- c(sprintf("%s/atac_pbmc_10k_nextgem_fragments.tsv.gz", data_dir)
              )
names(inputFiles) <- c(
    'pbmc_10k_atac'
    )

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 10, # Be careful about the minTSS and minFrags parameters. Start with 1, and 500
  minFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
)


# Create project
proj_name <- "pbmc_10x_atac"
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

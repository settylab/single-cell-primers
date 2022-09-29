library(ArchR)
set.seed(1)
system2("ml", "MACS2/2.2.6-foss-2019b-Python-3.7.4")
# Configure
addArchRThreads(threads = 12)
addArchRGenome('hg38')


# ################################################################################################
# Arrow files and project 

# Input files
data_dir = '/fh/fast/setty_m/grp/lab-datasets/bonemarrow-tcell-dep-multiome/cr-arc-results'
output_dir = '/fh/fast/setty_m/user/cjordan2/repositories/single-cell-primers/tcell-multiome-data'
dir.create(sprintf("%s/ArchR", output_dir))
setwd(sprintf("%s/ArchR", output_dir))
inputFiles <- c(sprintf("%s/rep1/atac_fragments.tsv.gz", data_dir)
              , sprintf("%s/rep2/atac_fragments.tsv.gz", data_dir))
names(inputFiles) <- c('tcell_dep_multiome_1', 'tcell_dep_multiome_2'
                       )

# Create Arrow files
#ArrowFiles <- createArrowFiles(
#  inputFiles = inputFiles,
#  sampleNames = names(inputFiles),
#  minTSS = 1,
 # Dont set this too high because you can always increase later
#  minFrags =3000, 
#  addTileMat = TRUE,
#  addGeneScoreMat = FALSE,
#  excludeChr = c('chrM')
#)
#output_dir = sprintf("%s/ArchR", output_dir)
#ArrowFiles = c(sprintf("%s/NA.arrow", output_dir), sprintf("%s/tcell_dep_multiome_1.arrow", output_dir), sprintf("%s/tcell_dep_multiome_2.arrow", output_dir), sprintf("%s/tcell_dep_multiome.arrow", output_dir))
# Create project
#proj_name <- "temp"
#proj <- ArchRProject(
  #ArrowFiles = ArrowFiles, 
  #outputDirectory = proj_name,
  #copyArrows = FALSE,
#)


# Subset of cells determined in RNA
# Multiome
#barcode_dir <- '/fh/fast/setty_m/user/cjordan2/repositories/single-cell-primers/data'
#multiome_cells = read.csv(sprintf("%s/tcell-dep_multiome_cells.csv", barcode_dir), stringsAsFactors=FALSE)[,2]
#multiome_cells <- intersect(multiome_cells, getCellNames(proj))
proj = loadArchRProject("/fh/fast/setty_m/user/cjordan2/repositories/single-cell-primers/tcell-multiome-data/ArchR/tcell-dep_multiome_atac")
# Project 
#proj_name = "tcell-dep_multiome_atac"
#proj <- subsetArchRProject(proj, multiome_cells, proj_name, force = FALSE)





# ################################################################################################
# PReprocesing

# SVD, Clustering, UMAP
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                       name = "IterativeLSI", scaleDims=FALSE, force=TRUE)#, varFeatures=100000)
var_features <- proj@reducedDims[["IterativeLSI"]]$LSIFeatures

# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
#proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', blacklist=blacklist, force=TRUE)



# Peaks 
#proj <- addClusters(input = proj,reducedDims = "IterativeLSI", force = TRUE)
##removed deprecated argument maxFragmentLength
#proj <- addGroupCoverages(proj, force = TRUE)
#proj <- addReproduciblePeakSet(proj, force = FALSE)
# Counts
#proj <- addPeakMatrix(proj,ceiling=10^9, force = TRUE)
##removed deprecated argument maxFragmentLength
proj <- addImputeWeights(proj)
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


#imputed gene_scores
imputed.weights <-getImputeWeights(proj)
write.csv(imputed.weights, sprintf('%s/export/imputed_weights.csv', proj_name), quote=FALSE)
imputed_scores <- imputeMatrix(mat = assay(gene.scores), imputeWeights = imputed.weights)
write.csv(imputed_scores, sprintf('%s/export/imputed_scores.csv', proj_name), quote=FALSE)

scores <-assays(gene.scores)
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



library(ArchR)
library(parallel)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(GenomicRanges)
library(plyranges)

set.seed(1)

# Configure
addArchRThreads(threads = 1)
addArchRGenome('hg38')


# ################################################################################################
# Arrow files and project 

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

# Subset of cells determined in RNA
# Multiome
multiome_path =  sprintf("%s/data", output_dir)
valid_barcodes = list()
rep1 = read.csv(sprintf("%s/tcell-dep_atac_cells_rep1.csv", multiome_path), stringsAsFactors = FALSE)[,2]
valid_barcodes[['rep1']] = rep1
rep2 = read.csv(sprintf("%s/tcell-dep_atac_cells_rep2.csv", multiome_path), stringsAsFactors = FALSE)[,2]
valid_barcodes[['rep2']] = rep2
##list item names must be replicate names. Each CSV file needs a header line. barcodes should not contain any replicate names.
setwd(sprintf("%s/ArchR", output_dir))


db <- makeTxDbFromGFF(file = '/fh/fast/setty_m/grp/tools/cellranger/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz', format= 'gtf', organism = 'Homo sapiens')

annotation <- createGeneAnnotation(TxDb = db, OrgDb = org.Hs.eg.db)
annotation_genes_filtered = annotation$genes
annotation_genes_filtered <- annotation_genes_filtered %>%
    filter(seqnames %in% chrs)
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

# Create Arrow files 
# Note that the TSS and Frags filter might result in some cells not being included. 
# Set these to 0 if you would like all cells to included.

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0, 
  minFrags = 0, 
  validBarcodes = valid_barcodes,
  geneAnnotation = annotation,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
  force = FALSE
)

proj_name <- "tcell_dep_ncrna"
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = proj_name,
  copyArrows = FALSE,
)


n_cells = nCells(proj)
print(sprintf("Number of total cells in project: %s", n_cells))
proj <- saveArchRProject(ArchRProj = proj)

# ################################################################################################
# PReprocesing

# SVD, Clustering, UMAP

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                       name = "IterativeLSI", scaleDims=TRUE, force=FALSE, outlierQuantiles = NULL) #maybe remove outlier quantiles
var_features <- proj@reducedDims[["IterativeLSI"]]$LSIFeatures

# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)

zproj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', genes = annotation_genes_filtered, threads = 1, force = TRUE)


# Peaks 
proj <- addClusters(input = proj,reducedDims = "IterativeLSI", force = TRUE)
proj <- addGroupCoverages(proj, maxFragmentLength = 147, force = TRUE, threads = 1)
proj <- addReproduciblePeakSet(proj, force = TRUE)
# Counts
proj <- addPeakMatrix(proj,ceiling=10^9, maxFragmentLength = 147, force = TRUE)

#UMAP
proj <- addUMAP(proj)

# ChromVAR
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

# Motif ranking
var_deviations <- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE)
markerMotifs <- sprintf("z:%s", var_deviations$name)
warnings()

# Motif scores
motifScores <- ArchR:::.getMatrixValues(proj, name = markerMotifs, matrixName = 'MotifMatrix', log2Norm = FALSE)

# Save 
proj <- saveArchRProject(ArchRProj = proj)


# ################################################################################################
proj_name = sprintf("%s/ArchR/tcell_dep_ncrna", output_dir)
# Export
dir.create(sprintf("%s/export", proj_name))
write.csv(motifScores, sprintf('%s/export/motif_scores.csv', proj_name), quote=FALSE)


message("successfully ChromVAR exported")
write.csv(getReducedDims(proj), sprintf('%s/export/svd.csv', proj_name), quote=FALSE)
write.csv(getCellColData(proj), sprintf('%s/export/cell_metadata.csv', proj_name), quote=FALSE)


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


scores <-assays(gene.scores)
# chromVar
motifScores <-   getMatrixFromProject(proj, 'MotifMatrix')
scores <- as.data.frame(assays(motifScores)[['z']])
write.csv(scores, quote=FALSE,
         sprintf("%s/export/chromvar_motif_scores.csv", proj_name))
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

proj <- saveArchRProject(ArchRProj = proj)

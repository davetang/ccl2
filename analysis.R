#load CAGEr library
library(CAGEr)

#load full human genome (UCSC version hg19)
library(BSgenome.Hsapiens.UCSC.hg19)

#object for all bam files
bam_file <- list.files('.', "*F512.bam$")

#remove bad library
bam_file <- bam_file[-3]

#create CAGEr object
cage_bam <- new("CAGEset",
                  genomeName = "BSgenome.Hsapiens.UCSC.hg19",
                  inputFiles = bam_file,
                  inputFilesType = 'bam',
                  sampleLabels = c('hIPS_1', 'hIPS_2', 'hIPS_ccl2_1',  'hIPS_ccl2_2', 'hIPS_ccl2_3')
                )

#CAGE defined TSSes are called CTSS
#HelicosCAGE reads don't have quality scores
#sequencingQualityThreshold=-1 parameter circumvents the lack of quality scores
getCTSS(cage_bam, mappingQualityThreshold=10, sequencingQualityThreshold=-1)
cage_bam_ctss <- CTSStagCount(cage_bam)

#tpm normalisation
normalizeTagCount(cage_bam, method = "simpleTpm")

#load library for parallelisation
library(multicore)

#tag clustering with 6 cores
clusterCTSS(object = cage_bam,
            threshold = 1,
            thresholdIsTpm = TRUE,
            nrPassThreshold = 1,
            method = "distclu",
            maxDist = 20,
            removeSingletons = TRUE,
            keepSingletonsAbove = 5,
            useMulticore = T,
            nrCores = 6
           )

#cumulative sum of CAGE signal along genomic region
cumulativeCTSSdistribution(cage_bam, clusters = "tagClusters", useMulticore = T, nrCores = 6)
quantilePositions(cage_bam, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

#aggregating tag clusters across multiple CAGE datasets
aggregateTagClusters(cage_bam, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)

#Extracting consensus clusters from CAGEset object
cage_bam_consensus_cluster <- consensusClusters(cage_bam)

#CAGE data based expression clustering
getExpressionProfiles(cage_bam, what = "consensusClusters",
                      tpmThreshold = 10, nrPassThreshold = 1,
                      method = "som", xDim = 4, yDim = 2
                      )

#analysis for shifting promoter usage
cumulativeCTSSdistribution(cage_bam, clusters = "consensusClusters", useMulticore = T, nrCores = 8)
scoreShift(cage_bam, groupX = c('hIPS_1', 'hIPS_2'), groupY = c('hIPS_ccl2_1',  'hIPS_ccl2_2', 'hIPS_ccl2_3'), useMulticore = T, nrCores = 8)
cage_bam_shifting_promoter <- getShiftingPromoters(cage_bam, tpmThreshold = 5, scoreThreshold = 0.6)

#preparing expression table for edgeR
for_edger <- cage_bam

#no normalisation using CAGEr
#will use TMM using edgeR later
normalizeTagCount(for_edger, method = "none")
clusterCTSS(object = for_edger,
            threshold = 1, thresholdIsTpm = FALSE,
            nrPassThreshold = 1, method = "distclu",
            maxDist = 20, removeSingletons = FALSE,
            useMulticore = T, nrCores = 8)
aggregateTagClusters(for_edger, tpmThreshold = 0, qLow = NULL , qUp = NULL, maxDist = 100)
for_edger_count <- for_edger@consensusClustersTpmMatrix
for_edger_consensus_cluster <- consensusClusters(for_edger)
for_edger_consensus_cluster$id <- paste(for_edger_consensus_cluster$chr, for_edger_consensus_cluster$start, for_edger_consensus_cluster$end, for_edger_consensus_cluster$strand, sep="_")
d <- for_edger_count
rownames(d) <- for_edger_consensus_cluster$id
head(d)
                      sample
consensus.cluster      hIPS_1 hIPS_2 hIPS_ccl2_1 hIPS_ccl2_2 hIPS_ccl2_3
  chr1_10160_10161_+        0      0           0           1           0
  chr1_11854_11855_+        0      0           0           1           0
  chr1_20216_20217_+        1      0           0           0           0
  chr1_57934_57935_+        0      0           1           0           1
  chr1_88768_88769_+        0      0           0           0           1
  chr1_227833_227834_+      1      0           0           0           0

library(edgeR)
group <- c(rep("C",2),rep("T",3))
d <- DGEList(counts = d, group=group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d, verbose=T)
Disp = 0.09247 , BCV = 0.3041
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d)
summary(decideTestsDGE(de.tgw, p.value=0.01))
   [,1]
-1      34
0  3719817
1       77
topTags(de.tgw)
Comparison of groups:  T-C
                                logFC   logCPM        PValue           FDR
chr4_74301879_74301944_+    -6.556924 4.981371 2.751677e-107 1.023604e-100
chr2_168468160_168468205_+  -2.001768 6.932641  1.477002e-75  2.747170e-69
chr2_46524249_46525133_+     2.268422 6.693334  4.814475e-54  5.969834e-48
chr11_116700593_116700629_+ -8.678544 2.944637  6.850884e-35  6.371199e-29
chr13_96204918_96205208_+    1.723369 6.224111  1.645455e-32  1.128967e-26
chr14_94854932_94854984_-   -6.836681 2.907701  1.820950e-32  1.128967e-26
chrX_119149289_119149903_-   3.399951 4.886728  1.395220e-26  7.414454e-21
chr17_66595846_66597739_-    1.181361 6.370246  1.349262e-25  6.273947e-20
chr10_8096538_8096794_+      2.287188 5.384277  1.198765e-22  4.954798e-17
chr3_186330847_186330952_+  -6.222377 2.356282  1.642015e-22  6.108178e-17

save.image("analysis_leave_one_out.Robject")


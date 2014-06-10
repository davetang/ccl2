#Entire analysis performed on CCL2 CAGE libraries
#please email me@davetang.org to comment or ask questions
#save this script, analysis.R, in the same directory as the 6 aligned libraries
#which should be as BAM files
#the names of the BAM files need to end with *F512.bam
#F512 means reads not passing QC

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

#cluster CTSS using single lineage clustering
clusterCTSS(object = for_edger,
            threshold = 1, thresholdIsTpm = FALSE,
            nrPassThreshold = 1, method = "distclu",
            maxDist = 20, removeSingletons = FALSE,
            useMulticore = T, nrCores = 8)

#aggregarting tag clusters across multiple CAGE datasets
aggregateTagClusters(for_edger, tpmThreshold = 0, qLow = NULL , qUp = NULL, maxDist = 100)

#count table
for_edger_count <- for_edger@consensusClustersTpmMatrix

#creating identifier for each tag cluster
for_edger_consensus_cluster <- consensusClusters(for_edger)
for_edger_consensus_cluster$id <- paste(for_edger_consensus_cluster$chr,
					for_edger_consensus_cluster$start,
					for_edger_consensus_cluster$end,
					for_edger_consensus_cluster$strand,
					sep="_")

#set threshold for independent filtering
rowsum_threshold <- 20

#set p-value threshold for all analyses
cutoff <- 0.05

#preparing object for edgeR
d2 <- for_edger_count
rownames(d2) <- for_edger_consensus_cluster$id

#independent filtering
d2 <- d2[rowSums(d2)>rowsum_threshold,]

#load library
library(edgeR)
#2 vs 3 analysis
group <- c(rep("C",2),rep("T",3))
d2 <- DGEList(counts = d2, group=group)

#TMM normalisation
d2 <- calcNormFactors(d2)

#calculate dispersion
d2 <- estimateCommonDisp(d2, verbose=T)
d2 <- estimateTagwiseDisp(d2)

#differential expression test
de.tgw2 <- exactTest(d2)

#summary of differentially expressed tag clusters
summary(decideTestsDGE(de.tgw2, p.value=cutoff))
#   [,1]
#-1    272
#0  112928
#1     147

#create data frame for results
#store normalised counts
data2 <- d2$pseudo.counts

#store differential analysis results
data2 <- cbind(data2, de.tgw2$table)

#calculate FDR
data2$FDR <- p.adjust(data2$PValue, method='BH')

#store dispersion of each tag cluster
data2$tw_dis <- d2$tagwise.dispersion

#prepare coordinates of each tag cluster
data_coord2 <- matrix(data=unlist(strsplit(row.names(data2), split="_")),
                      nrow= length(row.names(data2)),
                      byrow=T)
data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)
names(data_coord2) <- c('chr','start','end','strand')

#store coordinates of each tag cluster
data2 <- cbind(data2, data_coord2)

#create column for differential expression status
#1 for DE and 0 for not
data2$de <- as.numeric(data2$FDR<cutoff)

#convert coordinates to numeric
data2$start <- as.numeric(data2$start)
data2$end <- as.numeric(data2$end)

#load libraries
library(GenomicRanges)
library(GenomicFeatures)

#create GRanges object for all tag clusters
data_grange2 <- with(data2,
                    GRanges(chr, IRanges(start,
                                         end,
                                         names=row.names(data2)
                                         ),
                            strand
                            )
                    )

#get Entrez gene annotation using biomaRt
library("biomaRt")
ensembl <- useMart("ensembl")
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

#obtain all Entrez genes on these chromosomes
my_chr <- c(1:22, 'M', 'X', 'Y')

#obtain Entrez genes
my_entrez_gene <- getBM(attributes='entrezgene',
                        filters = 'chromosome_name',
                        values = my_chr,
                        mart = ensembl)

#attributes to fetch
my_attribute <- c('entrezgene',
                  'chromosome_name',
                  'transcript_start',
                  'transcript_end',
                  'strand')

#fetch from biomaRt
my_entrez_gene_loci <- getBM(attributes=my_attribute,
                             filters = c('entrezgene', 'chromosome_name'),
                             values = list(entrezgene=my_entrez_gene$entrezgene, chromosome_name=my_chr),
                             mart = ensembl)

#convert strand information from -1/1 to -/+ respectively
my_entrez_gene_loci$strand <- gsub(pattern='-1',
                                   replacement='-',
                                   my_entrez_gene_loci$strand)
my_entrez_gene_loci$strand <- gsub(pattern='1',
                                   replacement='+',
                                   my_entrez_gene_loci$strand)

#region around gene to call a TSS
span <- 200

#create new data frame for TSSs
my_entrez_tss <- my_entrez_gene_loci

#positive strand
#adjust the end position first, because we need the start position in our calculations
my_entrez_tss[my_entrez_tss$strand=='+','transcript_end']   <- my_entrez_tss[my_entrez_tss$strand=='+','transcript_start']+span
my_entrez_tss[my_entrez_tss$strand=='+','transcript_start'] <- my_entrez_tss[my_entrez_tss$strand=='+','transcript_start']-span

#negative strand
my_entrez_tss[my_entrez_tss$strand=='-','transcript_start'] <- my_entrez_tss[my_entrez_tss$strand=='-','transcript_end']-span
my_entrez_tss[my_entrez_tss$strand=='-','transcript_end']   <- my_entrez_tss[my_entrez_tss$strand=='-','transcript_end']+span

#add a 'chr' into the chromosome_name
my_entrez_tss$chromosome_name <- gsub(pattern="^",
                                      replacement='chr',
                                      my_entrez_tss$chromosome_name)

#write more informative column names
names(my_entrez_tss) <- c('id','chr','start','end','strand')

#create a Granges object for the Entrez TSSs
entrez_grange <- with(my_entrez_tss,
                      GRanges(chr, IRanges(start,
                                           end,
                                           names=id
                                           ),
                              strand
                              )
                      )

#overlap Entrez TSSs with all tag clusters
tc_entrez2 <- findOverlaps(data_grange2, entrez_grange)

#store number of overlaps
tc_entrez_count2 <- countOverlaps(data_grange2, entrez_grange)

#store number of overlaps in the big table of results
data2$oc <- as.vector(tc_entrez_count2)

#create data frame that matches the Entrez TSS to the tag cluster
match_hit2 <- data.frame(names(data_grange2)[queryHits(tc_entrez2)],
                        names(entrez_grange)[subjectHits(tc_entrez2)],
                        stringsAsFactors=F
                        )

#name the columns
names(match_hit2) <- c('query','subject')

#remove duplicated entries
match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

#name the rows
row.names(match_hit2) <- match_hit2$query

#store Entrez gene in big table of results
data2$gene <- match_hit2[row.names(data2),2]

#obtain HUGO gene symbols for Entrez genes
entrez_symbol <- getBM(attributes=c('entrezgene','hgnc_symbol'),
		       filters = 'entrezgene',
		       values = match_hit2$subject,
		       mart = ensembl)

#remove blank entries
entrez_symbol <- entrez_symbol[-grep("^$", entrez_symbol$hgnc_symbol, perl=T),]

#name columns
names(entrez_symbol) <- c('subject', 'symbol')

#remove duplicated entries
entrez_symbol <- entrez_symbol[!duplicated(entrez_symbol$subject),]

#name rows
row.names(entrez_symbol) <- entrez_symbol$subject

#merge with big table of results
data2$gs <- entrez_symbol[data2$gene,2]

#Gene Ontology Enrichment Analysis
#load libraries
library(GO.db)
library("GOstats")

#tag clusters differentially expressed with the CCL2 condition
up_ccl2_2 <- row.names(subset(data2, FDR<cutoff & logFC>0 & oc>0))

#universal list
universe2 <- row.names(subset(data2, oc>0))

#tag clusters to Entrez gene
up_ccl2_entrez2 <- match_hit2[up_ccl2_2,2]
up_ccl2_entrez2 <- unique(up_ccl2_entrez2)
universe_entrez2 <- match_hit2[universe2,2]
universe_entrez2 <- unique(universe_entrez2)

#load library
library(org.Hs.eg.db)

#parameters for hypergeometric test
ccl2_bp_param <- new('GOHyperGParams',
              geneIds=up_ccl2_entrez2,
              universeGeneIds=universe_entrez2,
              ontology='BP',
              pvalueCutoff=cutoff,
              conditional=F,
              testDirection='over',
              annotation="org.Hs.eg.db"
             )

#perform the test
ccl2_bp <- hyperGTest(ccl2_bp_param)

#store results
ccl2_bp_result <- summary(ccl2_bp)

#FDR
ccl2_bp_result$adj_p_value <- p.adjust(ccl2_bp_result$Pvalue, method="BH", n=1910)

#repeat analysis but for bFGF libraries
up_bfgf <- row.names(subset(data2, FDR<cutoff & logFC<0 & oc>0))
up_bfgf_entrez <- match_hit2[up_bfgf,2]
up_bfgf_entrez <- unique(up_bfgf_entrez)

bfgf_bp_param <- new('GOHyperGParams',
              geneIds=up_bfgf_entrez,
              universeGeneIds=universe_entrez2,
              ontology='BP',
              pvalueCutoff=cutoff,
              conditional=F,
              testDirection='over',
              annotation="org.Hs.eg.db"
             )
bfgf_bp <- hyperGTest(bfgf_bp_param)
bfgf_bp_result <- summary(bfgf_bp)
bfgf_bp_result$adj_p_value <- p.adjust(bfgf_bp_result$Pvalue, method="BH", n=1191)

#For hypoxia related genes
#load libraries
library(org.Hs.eg.db)
library(GO.db)

#create gene ontology object
go_object <- as.list(org.Hs.egGO2EG)

#all GO terms related to hypoxia via AmiGO
#http://amigo.geneontology.org/cgi-bin/amigo/search.cgi?search_query=hypoxia&search_constraint=term&action=new-search
#find all Entrez genes associated with these GO IDs
a <- go_object['GO:0097411']
b <- go_object['GO:0001666']
c <- go_object['GO:0070483']
d <- go_object['GO:0071456']
e <- go_object['GO:1900037']
f <- go_object['GO:1900038']
g <- go_object['GO:1900039']
h <- go_object['GO:1990144']
i <- go_object['GO:1902071']
j <- go_object['GO:1902072']
k <- go_object['GO:1902073']
l <- go_object['GO:0061418']
m <- go_object['GO:0061428']
n <- go_object['GO:0061419']
o <- go_object['GO:2000777']

#store all
all <- list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)

#keep unique list of Entrez genes
all <- unique(unlist(all, use.names=F))

#store as data frame
all <- data.frame(gene=all)

#store only annoated tag clusters
data_annotated <- data2[!is.na(data2$gene),]

#results of hypoxia related genes
all_merge <- merge(x=all, y=data_annotated, by='gene')

#save analysis
save.image("analysis.Robject")

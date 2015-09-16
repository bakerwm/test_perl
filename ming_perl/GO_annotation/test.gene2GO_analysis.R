
###
setwd("/home/wangming/work/R_work/GO_analysis/wk_dirs/temp/GOdb_H37Rv")

## ----load modified functions of clusterProfiler
suppressPackageStartupMessages(require(clusterProfiler))
source("enrichMap-v2.R")
source("barplot-v2.R")
library(ggplot2)
library(reshape2)
library(igraph)
library(DOSE)

### gene table
gff <- "../NC_000962.gff"
### load genes
Gff2GeneTable(gff)
load("geneTable.rda")
head(geneTable, 4)
H37Rv_genes <- as.character(geneTable$GeneID)

### prepare gene table
da <- read.table("/home/wangming/work/R_work/GO_analysis/wk_dirs/temp/H37Ra_H37Rv_exp.diff", header = TRUE)
db <- da[, c(1, 5:6, 8:10)]
### log2 > 1, < -1
d2 <- db[(db$log2.fold_change. > 1 | db$log2.fold_change. < -1), ] # log2 >1, <-1
d4 <- db[(db$log2.fold_change. > 2 | db$log2.fold_change. < -2), ] # log2 >2, <-2

### find GeneID
df <- subset(geneTable, Locus %in% d2$test_id)
genes <- as.character(df$GeneID)
head(genes)

#############################
# For GO analysis
#############################

### build GOmap
godb <- "/home/wangming/work/R_work/GO_analysis/wk_dirs/temp/GOdb_H37Rv/NC_000962.GO.db"
### build
gomap <- read.table(godb, colClasses = "character",
                    col.names = c("entrezgene", "go_accession"))
gomap <- gomap[c("go_accession", "entrezgene")] # change the order of column in lastest version of clusterProfiler v2.2.3
dim(gomap)
head(gomap, 4)
buildGOmap(gomap)

### 1. GO classification
ggo1 <- groupGO(genes, organism = "MTB", ont = "MF", level = 4, readable = TRUE)
head(summary(ggo1), 4)

### 2. Enrichment analysis
### Hypergeometric model
ego1 <- enrichGO(genes, organism = "mtb", ont = "CC", pvalueCutoff = 1, 
                 qvalue = 1, readable = TRUE)
head(summary(ego1), 4)

### 
ego2 <- gseGO(genes, organism = "H37Rv", ont = "CC", nPerm = 10, minGSSize = 10,
              pvalueCutoff = 0.01, verbose = FALSE)
head(summary(ego2), 4)
#############################
### KEGG pathway enrichment
#############################
### perform online annotation for KEGG analysis
ge2 <- as.character(df$Locus)
kk1 <- enrichKEGG(ge2, organism = "mtu", pvalueCutoff = 1, qvalueCutoff = 1)
head(summary(kk1))
kk2 <- gseKEGG(ge2, organism = "mtu", nPerm = 1, minGSSize = 2, pvalueCutoff = 1, verbose = FALSE)
head(summary(kk2))

#############################
# output figures
#############################
### 1. barplot
barplot(ggo1, drop = TRUE, showCategory = 10)
barplot(ego1, drop = TRUE, showCategory = 10)

### 2. enrichMap
enrichMap(ego1, n = 20)

### 3. cnetplot
cnetplot(ego1, categorySize = "pvalue")
#cnetplot(ego1, categorySize = "pvalue", foldChange = geneList)
#cnetplot(kk, categorySize = "geneNum", foldChange = geneList)

### 4. gseaplot
#gseaplot(kk2, geneSetID = "hsa04145")

### pathview from pathview package
library(pathview)

mtu01100 <- pathview(gene.data = ge2, pathway.id = "mtu01100",
                     species = "mtu")


hsa04110 <- pathview(gene.data = geneList, pathway.id = "hsa04110",
                     species = "hsa", limit = list(gene = max(abs(geneList)), 
                                                   cpd = 1))

## ----YGC's example-------------------------------------------------------
spd <- c("SPD_0065", "SPD_0071", "SPD_0293", "SPD_0295", "SPD_0296", "SPD_0297", 
         "SPD_0327", "SPD_0426", "SPD_0427", "SPD_0428", "SPD_0559", "SPD_0560", 
         "SPD_0561", "SPD_0562", "SPD_0580", "SPD_0789", "SPD_1046", "SPD_1047", 
         "SPD_1048", "SPD_1050", "SPD_1051", "SPD_1052", "SPD_1053", "SPD_1057",
         "SPD_1326", "SPD_1432", "SPD_1534", "SPD_1582", "SPD_1612", "SPD_1613", 
         "SPD_1633", "SPD_1634", "SPD_1648", "SPD_1678", "SPD_1919")
kk1 <- enrichKEGG(spd, organism = "spd")
kk2 <- gseKEGG(spd, organism = "spd")
head(summary(kk1), 4)
head(summary(kk2))

###########################################################################
### ----clusterProfiler manual --------------------------------------------

### Convert IDs to GeneID (symbol, uniprot ID, others)
## ----bitr----------------------------------------------------------------
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2", 
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1", 
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1", 
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",  
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",  
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")
head(eg)

### User should provides an annotation package, both fromType and toType can accept any types that supported.

### User can use idType to list all supporting types.
idType("org.Hs.eg.db")

### one two multiple types
ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), annoDb="org.Hs.eg.db")
head(ids)

## ----groupGO-------------------------------------------------------------
library("DOSE")
data(geneList)
gene <- names(geneList)[abs(geneList) > 2]
head(gene)
ggo <- groupGO(gene     = gene,
               organism = "human",
               ont      = "BP",
               level    = 3,
               readable = TRUE)
head(summary(ggo))

## ----enrichGO------------------------------------------------------------
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                organism      = "human",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(summary(ego))

## ----gseGO---------------------------------------------------------------
ego2 <- gseGO(geneList     = geneList,
              organism     = "human",
              ont          = "CC",
              nPerm        = 100,
              minGSSize    = 120,
              pvalueCutoff = 0.01,
              verbose      = FALSE)
head(summary(ego2))

## ----enrichKEGG----------------------------------------------------------
kk <- enrichKEGG(gene         = gene,
                 organism     = "human",
                 pvalueCutoff = 0.05, 
                 readable     = TRUE)
head(summary(kk))

## ----gseKEGG-------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = "human",
               nPerm        = 100,
               minGSSize    = 120,
               pvalueCutoff = 0.01,
               verbose      = FALSE)
head(summary(kk2))

## ----enrichDAVID, eval=FALSE---------------------------------------------
##david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY")

### create figures
## ----barplot, fig.height=5, fig.width=6----------------------------------
barplot(ggo, drop=TRUE, showCategory=12)


## ----barplot-enrich, fig.height=5, fig.width=8---------------------------
barplot(ego, showCategory=8)


## ----enrichMap, fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, out.width="0.9\\textwidth", fig.pos="h"----
enrichMap(ego)


## ----cnetplot, fig.height=14, fig.width=14-------------------------------
cnetplot(ego, categorySize="pvalue", foldChange=geneList)


## ----cnetplot-KEGG, fig.height=14, fig.width=14--------------------------
cnetplot(kk, categorySize="geneNum", foldChange=geneList)


## ----gseaplot, fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8, out.width="0.6\\textwidth", fig.pos="h"----
gseaplot(kk2, geneSetID = "hsa04145")


## ----viewKEGG, eval=FALSE------------------------------------------------
## library("pathview")
## hsa04110 <- pathview(gene.data  = geneList,
##                      pathway.id = "hsa04110",
##                      species    = "hsa",
##                      limit      = list(gene=max(abs(geneList)), cpd=1))


## ----gcSample------------------------------------------------------------
data(gcSample)
lapply(gcSample, head)


## ----comparecluster------------------------------------------------------
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(summary(ck))


## ----formula-------------------------------------------------------------
## formula interface
mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
                            '100127206', '100128071'),
                   group = c('A', 'A', 'A', 'B', 'B', 'B'),
                   othergroup = c('good', 'good', 'bad', 'bad',
                                  'good', 'bad'))
xx.formula <- compareCluster(Entrez~group, data=mydf, fun='groupGO')
head(summary(xx.formula))

## formula interface with more than one grouping variable
xx.formula.twogroups <- compareCluster(Entrez~group+othergroup,
                                       data=mydf, fun='groupGO')
head(summary(xx.formula.twogroups))


## ----compareCluster, fig.height=8, fig.width=8---------------------------
plot(ck)

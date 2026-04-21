#################################
# R script for exploring gene expression in the single cell dataset.
#################################
# remove all the objects from the R session
rm(list=ls()) 
#set the working directory: 
# ***********update "WD" with the path to your working directory on your computer:
WD="PATH TO YOUR WORKING DIRECTORY"
setwd(WD)

#Download the R-object with raw single-cell RNA-Seq data to your working directory
#all data:
#https://cells.ucsc.edu/sea-anemone-atlas/Nv2/all/AllData.Robj"
#retractor muscle data:
#https://cells.ucsc.edu/sea-anemone-atlas/Nv2/retractor-muscle/
#pSC data: 
#https://cells.ucsc.edu/sea-anemone-atlas/Nv2/pSC/


#load the libraries that we use here | install these packages and dependencies first!
library("Seurat") # for working with the dataset
library("pals") # for the colour palettes
library("ggplot2") # for increased plotting functionality
library("readxl") # to use an excel spreadsheet to load a list of gene models of interest.
library("RColorBrewer")

#load the data. choose all data (Fig2A) or a subset (Fig 2B, Fig S4). 
load("AllData.Robj") 
load ("retractor-muscle.Robj")
load ("pSC.Robj")

#conserve an unmodified version of the dataset for security reasons
data1 <- AllData
#or
data1 <- `retractor muscle`
#or
data1 <- pSC

###
#Additional filtering only for pSC data:
#To generate the Dotplot for pSC in Figure S4, mitotic cells cells were dropped from the subset. 
data1<-SetIdent(data1,value='IDs')
data1<-subset(data1,idents=c("NPC.1", "NPC.2", "NPC.g"))
#drop empty clusters
data1$IDs<-droplevels(data1$IDs) 
#tabulate cells left
x=as.data.frame(table(data1$IDs))
data1<-SetIdent(data1,value='IDs')
#filter out any clusters with less than 10 cells
data1<-subset(data1,idents=levels(data1)[which(x$Freq>10)]) 
#drop empty clusters
data1$IDs<-droplevels(data1$IDs) 
data1<-SetIdent(data1,value = 'IDs')
levels(data1)
View(as.data.frame(table(data1$IDs)))
#End of additional filtering only for pSC data 
###

#set some nice colour palettes:
gene.cp=c('lightgrey',rev(brewer.pal(11 , "Spectral" ))) #rainbow gene expression palette
clust.cp.separate = unique (c(cols25(25),alphabet2(26),glasbey(32),alphabet(26))) #collated maximally unique colour palettes
clust.cp.graded = unique(c(stepped3(16),stepped(20),stepped2(20)))# graded colour palettes
LibCP = brewer.paired(21)
names(LibCP)=levels(SetIdent(data1,value='orig.ident'))


#Download the Excel file with gene names from Cole et al. (2024) Updated single cell reference atlas for the starlet anemone Nematostella vectensis. Front Zool 21, 8. 
#from https://static-content.springer.com/esm/art%3A10.1186%2Fs12983-024-00529-z/MediaObjects/12983_2024_529_MOESM3_ESM.xlsx
#Save the first sheet of the table as "NV2GenesFinZ.xlsx" in your working directory
genes<-readxl::read_xlsx('NV2GenesFinZ.xlsx')

#define genes of interest (goi) used in FIG2A and FIG2B
mygenes = c("NV2.21749","NV2.19965","NV2.2966","NV2.252","NV2.437","NV2.1449","NV2.8023","NV2.19847","NV2.6269","NV2.6489","NV2.2978","NV2.2977","NV2.10628","NV2.15308","NV2.11643","NV2.13851","NV2.13010","NV2.12106","NV2.24312","NV2.13963","NV2.14185","NV2.122","NV2.12016","NV2.7394","NV2.7393","NV2.12730","NV2.23832")
#goi used in FIG-S4
mygenes = c("NV2.21749","NV2.19847","NV2.6269","NV2.6489","NV2.2978","NV2.2977","NV2.10628","NV2.15308","NV2.11643","NV2.13851","NV2.13010","NV2.12106","NV2.24312","NV2.13963","NV2.14185","NV2.122","NV2.12016","NV2.7394","NV2.7393","NV2.12730","NV2.23832")

#index where these are
GOI = unique(genes$gene.name[match(mygenes,genes$geneID)])

# -----------------------------
# Step 1: Calculate average expression per group
# -----------------------------
avg <- AggregateExpression(
  data1,
  features = unique(GOI),
  group.by = "ID.separate" #for rm and pSC, group.by="IDs"
)$RNA

# Convert to matrix (important for downstream functions)
avg <- as.matrix(avg)

# -----------------------------
# Step 2: Optionally scale genes
# (clusters by expression pattern instead of magnitude)
# -----------------------------
avg_scaled <- t(scale(t(avg)))

gene_order <- unique(GOI) # or gene_order <- rev(unique(GOI))

# -----------------------------
# Step 3: Plot DotPlot with clustered gene order
# -----------------------------
library(dplyr)
library(tidyr)
library(readr)

dotPlotFigure <- DotPlot(
  data1,
  assay = "RNA",
  features = gene_order,
  scale.by = "size",
  scale.max = 100,
  col.min = 0,
  col.max = 3,
  group.by = "ID.separate", #for rm and pSC, group.by="IDs"
  cols = c("lightgrey","red")
) +
  RotatedAxis() +
  coord_flip() +
  theme(legend.position = "bottom")

renameTab <- read_tsv("geneID_newName.map")
new_order <- renameTab$new_gene_name[match(gene_order, renameTab$gene_short_name)]

dotPlotFigure$data <- dotPlotFigure$data %>% 
  left_join(renameTab, by=c("features.plot"="gene_short_name")) %>% 
  mutate(features.plot=new_gene_name) %>% 
  select(avg.exp,pct.exp,features.plot,id,avg.exp.scaled) %>%
  mutate(features.plot=factor(features.plot, levels=new_order))
print(dotPlotFigure)

# -----------------------------
# Step 4: All the raw expression values per cell are in data1, but you can export them as follows
# -----------------------------

library(dplyr)
library(tidyr)

expr <- FetchData(data1, vars = c("ID.separate", #for rm and pSC subset "IDs"
                                  "orig.ident", gene_order))

expr$cell_id <- rownames(expr)

expr_long <- expr %>%
  pivot_longer(
    cols = all_of(gene_order),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  rename(
    group = ID.separate, #for rm and pSC subset, group.by="IDs"
    sample_id = orig.ident
  ) %>%
  select(cell_id, sample_id, group, gene, expression)
#Export the values. Replace the FILENAME name below with an actual file name, 
#e.g. "Fig2A_Plot_raw_expression.csv"
write.csv(expr_long, "FILENAME_raw_expression.csv", row.names = FALSE)


# -----------------------------
# Step 5: Get summary statistics for what the plot shows
# ----------------------------- 
summary_stats <- expr_long %>%
  group_by(group, gene) %>%
  summarise(
    n_cells = n(),
    mean_expression = mean(expression),
    pct_expressing = mean(expression > 0) * 100,
    sd_expression = sd(expression),
    sem_expression = sd(expression) / sqrt(n_cells),
    .groups = "drop"
  )
#Export the values. Replace the FILENAME name below with an actual file name, 
#e.g. "Fig2A_Plot_summary_stats.csv"
write.csv(summary_stats, "FILENAME_summary_stats.csv", row.names = FALSE)
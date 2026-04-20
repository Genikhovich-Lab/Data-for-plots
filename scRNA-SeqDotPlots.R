#################################
# Here is an R script for exploring gene expression in the single cell dataset.
#################################
# remove all the objects from the R session
rm(list=ls()) 
#set the working directory: 
# ***********update "WD" with the path to your working directory on your computer:
WD="PATH TO YOUR WORKING DIRECTORY"
setwd(WD)

#Download the R-object with all the raw single-cell RNA-Seq data from
#https://cells.ucsc.edu/sea-anemone-atlas/Nv2/all/AllData.Robj"
#to your working directory


#load the libraries that we use here | install these packages and dependencies first!
library("Seurat") # for working with the dataset
library("pals") # for the colour palettes
library("ggplot2") # for increased plotting functionality
library("readxl") # to use an excel spreadsheet to load a list of gene models of interest.
library("RColorBrewer")

#load the data. 
load("AllData.Robj")

#conserve an unmodified version of the dataset for security reasons
data1 <- AllData 

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

#Download Excel files with a column geneID containing the NV2 gene model names 
#from https://github.com/Genikhovich-Lab/Data-for-plots into your working directory.
#Replace the FILENAME name below with an actual file name, e.g. "Plot1_IDs.xlsx":
excelgenes = read_excel('FILENAME.xlsx')
mygenes = excelgenes$geneID #ignore the rest of the worksheet
#index where these are
GOI = unique(genes$gene.name[match(mygenes,genes$geneID)])

# -----------------------------
# Step 1: Calculate average expression per group
# -----------------------------
avg <- AggregateExpression(
  data1,
  features = unique(GOI),
  group.by = "ID.separate"
)$RNA

# Convert to matrix (important for downstream functions)
avg <- as.matrix(avg)

# -----------------------------
# Step 2: Optionally scale genes
# (clusters by expression pattern instead of magnitude)
# -----------------------------
avg_scaled <- t(scale(t(avg)))

# -----------------------------
# Step 3: Cluster genes by similarity
# -----------------------------
gene_dist <- as.dist(1 - cor(t(avg_scaled)))

gene_clust <- hclust(gene_dist)

gene_order <- rownames(avg_scaled)[gene_clust$order]

# -----------------------------
# Step 4: Plot DotPlot with clustered gene order
# -----------------------------
DotPlot(
  data1,
  assay = "RNA",
  features = gene_order,
  scale.by = "size",
  col.min = 0,
  col.max = 3,
  group.by = "ID.separate",
  cols = c("lightgrey","red")
) +
  RotatedAxis() +
  coord_flip() +
  theme(legend.position = "bottom")


# -----------------------------
# Step 5: All the raw expression values per cell are in data1, but you can export them as follows
# -----------------------------

library(dplyr)
library(tidyr)

expr <- FetchData(data1, vars = c("ID.separate", "orig.ident", gene_order))

expr$cell_id <- rownames(expr)

expr_long <- expr %>%
  pivot_longer(
    cols = all_of(gene_order),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  rename(
    group = ID.separate,
    sample_id = orig.ident
  ) %>%
  select(cell_id, sample_id, group, gene, expression)
#Export the values. Replace the FILENAME name below with an actual file name, 
#e.g. "Plot1_raw_expression.csv"
write.csv(expr_long, "FILENAME_raw_expression.csv", row.names = FALSE)


# -----------------------------
# Step 6: Get summary statistics for what the plot shows
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
#e.g. "Plot1_summary_stats.csv"
write.csv(summary_stats, "FILENAME_summary_stats.csv", row.names = FALSE)

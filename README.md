# Data-for-plots
scRNA-SeqDotPlots.R contains the script to generate the expression dot plots and, if necessary, extract raw data and expression stats from publicly accessible scRNA-Seq datastet published in https://doi.org/10.1186/s12983-024-00529-z. All sequencing data from this paper is available as an R object at https://cells.ucsc.edu/sea-anemone-atlas/Nv2/all/AllData.Robj
Please note that downloading it directly into R-Studio failed in our hands. Save it locally and proceed as described in the script. 
The .xlsx contain geneIDs for the corresponding plots.
The geneID_newName.map file is a tab-delimited text file allowing to append NV2 gene model numbers to the gene names from AllData.Robj

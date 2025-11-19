setwd("D:/R/0902/monocle")
library(monocle)

rm(list = ls())

matrix <- Treg.intergrated@assays$RNA@data
matrix <- as.matrix(matrix)
pd <- new("AnnotatedDataFrame",data = Treg.intergrated@meta.data)
fdata <- data.frame(gene_short_name = rownames(matrix),row.names = row.names(matrix))
fd <- new("AnnotatedDataFrame",data = fdata)
monocle <- newCellDataSet(matrix,phenoData = pd,featureData = fd,
                          lowerDetectionLimit = 0.2,
                          expressionFamily = negbinomial.size())
monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)
monocle <- detectGenes(monocle, min_expr = 0.05)
expressed_genes <- row.names(subset(fData(monocle),
                                    num_cells_expressed >= 500))
print(head(pData(monocle)))
diff_test_res <- differentialGeneTest(monocle[expressed_genes,],fullModelFormulaStr="~seurat_clusters")
write.csv(diff_test_res,file="diff_test_res.csv")

#
monocle_ordering_genes <-
  row.names(diff_test_res)[order(diff_test_res$qval)][0:500]

#
monocle <- setOrderingFilter(monocle,
                             ordering_genes = monocle_ordering_genes)
#
monocle <- reduceDimension(monocle, method = 'DDRTree',
                           max_components = 2,
                           num_dim = 6,
                           reduction_formula = "~ orig.ident")
#
monocle <- orderCells(monocle)

plot_cell_trajectory(monocle, cell_size = 3 , color_by = "State")

ggsave("monocle.png",width = 7 , height = 5, dpi=600)

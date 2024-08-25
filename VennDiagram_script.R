# Draw VennDiagram and Extract common genes
#install.packages("VennDiagram")
library(VennDiagram)
??VennDiagram   # help phage

# Upload data
dataset1 = read.delim(file = "MCC.txt", header = T, sep = '\t')
SET1 = dataset1$Gene.symbol   # Keep gene symbol in a variable

dataset2 = read.delim(file = 'MC.txt', header = TRUE, sep = '\t')

SET2 = dataset2$Gene.symbol


v1 <- venn.diagram(list(GSE39612=SET1, GSE15605=SET2),
                   fill = c("red", "yellow"),
                   alpha = c(0.5, 0.5),
                   filename=NULL,
                   na= 'stop') 

grid.draw(v1)       # Draw figure


# Customize our Venn Diagram

v2 <- venn.diagram(list(GSE39612=SET1, GSE15605=SET2),
                   fill = c("green", "yellow"),
                   alpha = c(0.5, 0.5),
                   filename=NULL,
                   na= 'stop',
                   height = 2000, 
                   width = 2000, resolution = 500, imagetype = "tiff",
                   main = 'This is venn diagram',
                   sub = "Common gene") 
grid.newpage()        # Open new page
grid.draw(v2)       # draw figure





v3 <- venn.diagram(list(SET1, SET2),
                   fill = c("red", "green"),
                   category.names = c('GSE39612', 'GSE15605'),
                   cat.cex= c(1,1),
                   cat.col= c('red','green'),
                   alpha = c(0.5, 0.5),
                   filename=NULL,
                   na= 'stop',
                   height = 2000, 
                   width = 2000, 
                   resolution = 500, 
                   imagetype = "tiff",
                   main = 'This is venn diagram',
                   sub = "Common gene") 


grid.newpage()        # Open new page
grid.draw(v3)       # draw figure


# Extract Common gene by using intersect function.

common_gene= intersect(SET1, SET2)


Common_gene=as.data.frame.factor(common_gene)
Common_gene

# or

comgene = dataset1[dataset1$Gene.symbol %in% common_gene, c("Gene.symbol", "Gene.title")]  
write.csv(comgene, file = 'Commongene.csv')


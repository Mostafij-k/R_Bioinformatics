# Draw heatmap

#install.packages("RColorBrewer")
#install.packages("pheatmap")

library(ggplot2) 		
library(RColorBrewer)  # Help to extract different col0r
library(pheatmap)

# Upload Data
dataset= read.csv('heatmap_adj.csv', header = TRUE, sep = ',')
row.names(dataset)= dataset$Gene
colnames(dataset)= c('Gene', 'COVID-19', 'MCC','MC')
dataset= dataset[, -1]
dataset= data.matrix(dataset)   # Convert dataframe to matrix

        # Construct heatmap
help(pheatmap)  # open help page

pheatmap(dataset, border_color = 'red',
         cellwidth = 40,
         fontsize_col = 15,
         scale='column')

C= colorRampPalette(rev(brewer.pal(n = 7, name =
                                              "RdYlBu")))(100)  # Default color for pheatmap

pheatmap(dataset, border_color = C,
         cellwidth = 40,
         fontsize_col = 15,
         scale='column')



# Use different color
display.brewer.all()   # It shows all available color palette 


greys= colorRampPalette(brewer.pal(9, 'Greys'))(100)

pheatmap(dataset, border_color = 'red',
         cellwidth = 40,
         fontsize_col = 15,
         scale='column',
         color = greys)

pairs= colorRampPalette(brewer.pal(9, 'Paired'))(100)

pheatmap(dataset, border_color = 'red',
         cellwidth = 40,
         fontsize_col = 15,
         scale='column',
         color = pairs)

purd= colorRampPalette(brewer.pal(7, 'PuRd'))(100)


pheatmap(dataset, border_color = 'red',
         cellwidth = 40,
         fontsize_col = 15,
         scale='column',
         color = purd)


blue= colorRampPalette(brewer.pal(7, 'Blues'))(100)


pheatmap(dataset, border_color = 'red',
         cellwidth = 40,
         fontsize_col = 15,
         scale='column',
         color = blue)



pheatmap(dataset, border_color = 'red',
         cellwidth = 40,
         fontsize_col = 18,
         color = pairs,
         scale='column',
         cluster_rows = FALSE,
         display_numbers = TRUE)









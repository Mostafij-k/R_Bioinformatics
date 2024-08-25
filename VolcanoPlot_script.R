# Draw Volcano plot
# Install your requre packages 

#install.packages("ggplot2")
#install.packages("plotly")

library(ggplot2)
library(plotly)

dataset=read.delim(file = 'MC.txt', header = TRUE, sep='\t')

dataset1=as.data.frame(dataset)  # Functions to check if an object is a data frame, or coerce it if possible

View(dataset1)



               # Volcano plot


ggplot(dataset1, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point() #simple way of volcano plot(){Basic}

# P value and adj p value both are indicator of statistical significance
# lower p value is more significance
# .05 or 5% is threshold of p value
# lower .05 means more significance and over .05 lower significance

# Add a column of NAs
dataset1$Category <- "Not Sig"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 

dataset1$Category[dataset1$logFC > 2 & dataset1$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.05, set as "DOWN"

dataset1$Category[dataset1$logFC < -2 & dataset1$adj.P.Val < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"

p <- ggplot(data=dataset1, aes(x=logFC, y=-log10(adj.P.Val), col=Category)) +
  geom_point(size = 1.5, alpha =1 , na.rm = T, shape = 21)   # na.rm = NA remove



# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-2,2), col="red",linetype="dashed") +
  theme_bw(base_size = 25) +         # theme_bw will create border (fixed size)
  theme(legend.position = "right") +
  labs(x="log2 fold change",
       y="-Log10 adj p value",
       title="GSE15605(MC)")+
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed")+
  xlim(-6,6)+
  ylim(0,10)   



ggplotly(p2)

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
#p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))  # Create your own discrete scale

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "Not Sig")
p3 <- p2 + scale_colour_manual(values = mycolors)
ggplotly(p3)






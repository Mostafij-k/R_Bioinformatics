# Microarray data analysis Using R

         # We will analyze Markel cell carcinoma (GSE39612)

# Install require packages

# BiocManager::install("Biobase")
# 
# BiocManager::install('GEOquery')
# 
# BiocManager::install("limma")
# 
# install.packages('ggplot2')
# 
# install.packages("plotly")


library(Biobase) # Biobase contains standardized data structures to represent genomic data.

library(GEOquery)   # Download data set and process it

library(limma)      # Differential gene expression analysis

library(ggplot2)

library(plotly)

# Load series and platform data from GEO

gset <- getGEO("GSE39612", GSEMatrix =TRUE, AnnotGPL=TRUE)   
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1  
gset <- gset[[idx]]  # It convert large list to large expression Set

# getGEO function use to download relevant data set
# If data set contain more platform, you will be asked to select desire platform
# # Platform basically chip type
# Large list contains different type of information
# Expression set = It is designed to combine several different sources of information into a single convenient structure.
# S4 is more complicated object in R that is used by various bioconductor packages.
# S4 object helps to store lot of information in single object.

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
make.names(fvarLabels(gset))
# gset contains column names. This function includes this names in our final output table called toptable


# group membership for all samples
gsms <- "XXXXXX111111111111111111111111111111111100000000000000000000000000000000000000000000000000000000000000001111111111111111111111111111111111"
sml <- strsplit(gsms, split="")[[1]]  # It split GSM and take one by one

# filter out excluded samples (marked as "X")
sel <- which(sml != "X") # this function excludes X from gsms
sml <- sml[sel] 
gset <- gset[ ,sel]  # subset
# So now gset contain only our desire sample

# log2 transformation
ex <- exprs(gset)  # This function give expression value for sample of interest
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)) # various quantile
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) #  # This is the of checking if the data have already been log transformed
LogC
# Log transformation is a data transformation method in which it replaces each variable x with a log(x).
# na.rm = tells the function whether or not to remove NA values from the calculation.
# Limma package deserve log transform data
# In simple terms, a quantile is where a sample is divided into equal-sized = meadian

if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# Assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("control","Case"))  # Make valid names out of character vectors.
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)  # this code is a way to convert a factor to a matrix
colnames(design) <- levels(gs) # set the names to columns of a matrix.

# design is a mathematical expression of our experiment (define experiment)
# gset ~ group = gset depend on group
# + 0 means that the model will not have an intercept.
# Factor in R is a variable used to categorize and store the data. Factor in R is also known as a categorical variable


fit <- lmFit(gset, design)  # fit linear model
# linear model makes relationship between dependent variable and independent variable

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")  # Takes multiple elements from the multiple vectors and concatenates them into a single element.
cont.matrix <- makeContrasts(contrasts=cts, levels=design)  # It specifies what comparison, we want to make
fit2 <- contrasts.fit(fit, cont.matrix)   # Compute Contrasts from Linear Model Fit

# A negative coefficient suggests that as the independent variable increases, the dependent variable tends to decrease.



# compute statistics and table of top significant genes
fit2 <- eBayes(fit2) # these functions are used to rank genes in order of evidence for differential expression.
data <- topTable(fit2, adjust="fdr", sort.by="B", number=60000)  # Extract a table of the top-ranked genes from a linear model fit.


data2 <- subset(data, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

## Save our data with logFC > 2 and P.value 

mydata=data2[data2$P.Value<.05,]        
mydata2=mydata[abs(mydata$logFC)>2,]
mydata2[mydata2==''] = NA
mydata3=na.omit(mydata2)


write.table(mydata3, file='MCC2.txt', row.names=F, sep="\t")



ggplot(data2, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point() #simple way of volcano plot(){Basic}



# Add a column of NAs
data2$Category <- "Not Sig"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 

data2$Category[data2$logFC > 2 & data2$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.05, set as "DOWN"

data2$Category[data2$logFC < -2 & data2$adj.P.Val < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"

p <- ggplot(data=data2, aes(x=logFC, y=-log10(adj.P.Val), col=Category)) +
  geom_point(size = 1.5, alpha =1 , na.rm = T, shape = 21)   # na.rm = NA remove



# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-2,2), col="red",linetype="dashed") +
  theme_bw(base_size = 25) +         # theme_bw will create border (fixed size)
  theme(legend.position = "right") +
  labs(x="log2 fold change",
       y="-Log10 adj p value",
       title="GSE39612(MCC)")+
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed")+
  xlim(-6,6)+
  ylim(0,100)   



ggplotly(p2)

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
#p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))  # Create your own discrete scale

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "Not Sig")
p3 <- p2 + scale_colour_manual(values = mycolors)
ggplotly(p3)


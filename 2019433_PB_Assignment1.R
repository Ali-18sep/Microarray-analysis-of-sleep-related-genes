library(GEOquery) # importing libraries
library(limma)
library(umap)

library(DataExplorer)

#downloading the Dataset
gset <- getGEO("GSE58309", GSEMatrix=TRUE, AnnotGPL=TRUE)
if(length(gset)>1) idx<-grep("GPL6244", attr(gset, "names")) else idx<-1
gset <- gset[[idx]]


fvarLabels(gset) <- make.names(fvarLabels(gset))
gsms <- paste0("0101010101")
sml <- strsplit(gsms, split="")[[1]]

sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


#log Transformation
ex <- exprs(gset)
qx <- as.numeric((quantile(ex, c(0,0.25,0.5,0.75, 0.99, 1.0), na.rm = T)))

logC <- (qx[5]>100) || (qx[6]-qx[1]>50 && qx[2]>0)
if(logC) {ex[which(ex<=0)]<- NaN
ex<- log2(ex)}

logC <- (qx[5]>100) || (qx[6]-qx[1]>50 && qx[2]>0)
if(logC) {ex[which(ex<=0)]<- NaN
exprs(gset) <- log2(ex)}


# box & whisker plot
par(mar=c(7,4,2,1))
title <- paste ("GSE58309", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE58309", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# Histogram plot

plot_histogram(ex)


gs <- factor(sml)
groups <- make.names(c("S1","F1"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model


cts <- paste0(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)


fit2 <- eBayes(fit2, 0.01)
num_of_row <- nrow(ex)
print(num_of_row)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
DE <- tT[tT$adj.P.Val<0.05,]


#Writing the Gene_symbol of differentiated gene in the file R_Assignment1.csv
write.table( DE$Gene.symbol, file = "R_Assignment1.csv", sep = "\t" , row.names = FALSE, col.names = FALSE)


#R_Assignment1



#exPlt <- rnorm(ex,n=1000, m=80, sd=20)

#GPL1261

#install.packages("DataExplorer")


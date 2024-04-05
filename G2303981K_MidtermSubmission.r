# Discrete distributions
#
# Given that sequencing technology is very accurate (0.9999), compute the probability of seeing 10 errors in a sequence of 1000 bp.
#
# The following function is used to determine the probably of observing a certain number of errors (number_of_errors) in a sequence of a certain length (length_of_seq) with an error rate (error_rate). Fill up the function to return the correct probability.

# function to get the error probability
getErrorProbability <- function(number_of_errors, length_of_seq, error_rate=0.9999){
  # TODO: Add code to this function to return the error probability (Q1 3 marks)
  prob<-(choose(length_of_seq,number_of_errors))*((1-error_rate)^number_of_errors)*(error_rate^(length_of_seq-number_of_errors))
  return(prob)
}
P<-getErrorProbability(10,1000,error_rate=0.9999)

 # write out the answer to a text file
write(format(getErrorProbability(10, 10000, error_rate=0.9999), digits=3), file="q1.txt")

# Visualization
#
# The following code generates a test data.frame d and uses ggplot2 to plot a box plot.
#
# The data.frame contains of a column x with 3 values A, B and C. Column y contains numerical data. Create a box plot with the x in the X axis with y as the Y axis. There should be 3 individual box plots, one for each unique value of x (A, B, C). Box plot A should be filled in red, box plot B should be filled in blue and box plot C should be filled in green. The plot should have a title with the text "Box plot".

library(ggplot2)

set.seed(1027820)
d <- data.frame(x=sample(c("A", "B", "C"), size=1000, replace=T), y=rnorm(1000, mean=100, sd=10))

pdf("q2.pdf", width=8, height=6)
# TODO: complete the following line of code to generate a box plot. (Q2 3 marks)
p <- ggplot(d, aes(x, y)) + 
  geom_boxplot(fill=c("red","blue","green"))
print(p)
dev.off()

# PCA
#
# The following code generates a matrix of gene expression values with samples as the columns and genes as the rows. There are 30 samples in 3 sample groups prefixed with A, B and C

set.seed(19538282)
d <- rbind(
  t(sapply(rnorm(900, mean=7, sd=1.5), function(x) rnorm(30, mean=x, sd=0.25))),
  cbind(
    t(sapply(rnorm(100, mean=7, sd=1.5), function(x) rnorm(10, mean=x, sd=0.25))),
    t(sapply(rnorm(100, mean=7, sd=1.5), function(x) rnorm(10, mean=x, sd=0.25))),
    t(sapply(rnorm(100, mean=7, sd=1.5), function(x) rnorm(10, mean=x, sd=0.25)))
  )
)
colnames(d) <- as.character(sapply(c("A", "B", "C"), function(x) paste0(x, 1:10)))
rownames(d) <- paste0("Gene", 1:nrow(d))

# Complete the following function to generate and return the prcomp return value for a given matrix with genes as the rows and samples as the columns. The PCA should perform the dimension reduction on the genes in order to cluster the samples.

# function to return the PCA object
getPCAObject <- function(input_matrix) {
  # TODO: complete the code to return the prcomp return object. (Q3 3 marks)
  return(prcomp(input_matrix))
}

# run the PCA
fit <- getPCAObject(d)

# Plot the PCA with PC1 on the axis and PC2 on the Y axis, the X and Y axis should be labelled with the PC and the variance accounted for.

pdf("q4.pdf", width=8, height=6)
# TODO: complete the following line of code to generate a PCA plot with PC1 on the axis and PC2 on the Y axis. The X and Y axis should be labelled with the PC and the variance accounted for. There should be a title with the text "PCA". The sample groups A, B and C should be colored differently. (Q4 4 marks)
PC1<-fit$x[,1]
PC2<-fit$x[,2]
p <- ggplot(,aes(x=PC1, y=PC2, label=format(fit$sdev[1],fit$sdev[2],digits=3))) + geom_point(size=3) + labs(title="PCA")
print(p)
dev.off()

# Fill up the codes in the following 3 functions to generate 3 text files, one for the PC, one for the variance accounted for and one for the loadings.

getPC <- function(pca_obj) {
  SampleID <- data.frame(rownames(pca_obj$x))
  data <- cbind(SampleID,as.data.frame(pca_obj$x,row.names = F))
  colnames(data)[1]<-"Sample"
  return(data)
  # TODO: write codes to return a data.frame with the first column being the sample ID followed by the PC (from PC1 to PC30). The sample column should be entitled "Sample" and the PC columns should be in the format of PC1. Rows should be ordered from A1 to C10. PC columns should be ordered from PC1 to PC30. (Q5 2 marks)
}
write.table(getPC(fit), file="q5.tsv", quote=F, na="", sep="\t", col.names=T, row.names=F)

getVariance <- function(pca_obj) {
  Variance<-data.frame(pca_obj$sdev)
  PC<-data.frame(colnames(pca_obj$x))
  data<- cbind(PC,Variance)
  colnames(data)<-c("PC","Variance")
  return(data)
  # TODO: write codes to return a data.frame with the first column being the PC (PC1 to PC30) and the second column being the variance accounted for in percentage. The column names should be PC and Variance. The PC column should be ordered from PC1 to PC30. (Q6 2 marks)
}
write.table(getVariance(fit), file="q6.tsv", quote=F, na="", sep="\t", col.names=T, row.names=F)

getLoadings <- function(pca_obj) {
  Gene<-data.frame(rownames(pca_obj$x))
  PC<-data.frame(pca_obj$x)
  PC<-PC[,order(as.numeric(substring(colnames(PC),3,nchar(colnames(PC)))))]
  data<- cbind(Gene,PC)
  data<-data[order(as.numeric(substring(data[,1],5,nchar(data[,1])))),]
  return(data)
  # TODO: write codes to return a data.frame with the first column being the Gene (Gene1 to Gene1000) followed by the PC (from PC1 to PC30). The values for the PC columns are the loadings for each gene to each PC. The PC columns should be ordered from PC1 to PC30. The genes should be ordered from Gene1 to Gene1000. (Q7 2 marks)
}
write.table(getLoadings(fit), file="q7.tsv", quote=F, na="", sep="\t", col.names=T, row.names=F)

# Clustering

# The following code generates a matrix of gene expression values with samples as the columns and genes as the rows. There are 30 samples in 3 sample groups prefixed with A, B and C

set.seed(19538282)
d <- rbind(
  t(sapply(rnorm(900, mean=7, sd=1.5), function(x) rnorm(30, mean=x, sd=0.25))),
  cbind(
    t(sapply(rnorm(100, mean=7, sd=1.5), function(x) rnorm(10, mean=x, sd=0.25))),
    t(sapply(rnorm(100, mean=7, sd=1.5), function(x) rnorm(10, mean=x, sd=0.25))),
    t(sapply(rnorm(100, mean=7, sd=1.5), function(x) rnorm(10, mean=x, sd=0.25)))
  )
)
colnames(d) <- as.character(sapply(c("A", "B", "C"), function(x) paste0(x, 1:10)))
rownames(d) <- paste0("Gene", 1:nrow(d))

# Complete the following function to perform a hierarchical clustering on the generated data and return the hclust return value for a given matrix with genes as the rows and samples as the columns. The function should cluster the samples.

getHclust <- function(input_matrix) {
  data<-hclust(dist(input_matrix, method = "euclidean", diag = FALSE, upper = FALSE, p = 2), method = "aver")
  return(data)
  # TODO: write the codes to run a hierarchical clustering using Euclidean distances and average linkage. (Q8 3 marks)
}
fit <- getHclust(d)

pdf("q9.pdf", width=8, height=6)
library(ggplot2)
#BiocManager::install("ggdendro")
library(ggdendro)
ggdendrogram(fit, rotate = FALSE, size = 2)
# TODO: complete the following line of code to plot the dendrogram of the hierarchical clustering (Q9 2 marks)
plot()
dev.off()

# Complete the code in the following function to run a kmeans clustering with 3 centers on the generated data. Use the "kmeans" function in R with iter.max=10000. The function should return a data.frame with 2 columns, the first is the Sample and the second is the Cluster number (1-3). The column names of the data.frame should be "Sample" and "Cluster". The samples should be sorted from A1 to C10.
getKMeans <- function(input_matrix) {
  k<-data.frame(input_matrix)
  Gene<-data.frame(rownames(k))
  data<-kmeans(k, centers = 3, nstart = 20)
  Cluster<-as.data.frame(data$cluster)
  data<-cbind(Gene,Cluster)
  colnames(data)<-c("Gene","Cluster")
  return(data)
  # TODO: write codes to perform a kmeans clustering with 3 centers. (Q10 3 marks)
}
set.seed(1925623)
write.table(getKMeans(d), file="q10.tsv", quote=F, na="", sep="\t", col.names=T, row.names=F)

# Hypothesis testing

# Complete the code in the following function to return True if the data is normally distributed and False otherwise using the Shapiro-Wilk test

isNormal <- function(values) {
  sp<-shapiro.test(values)
  if (sp$p.value>0.05) {print("TRUE")} 
    else {print("FALSE")}
  # TODO: write codes to return True when the data is normally distributed and False otherwise. (Q11 2 marks)
}
write(isNormal(as.numeric(d[1, ])), file="q11.txt")

# The following codes generate a data.frame of gene expression data.

set.seed(195382)
d <- cbind(
  sample=paste0("S", 1:30),
  sample_group=as.character(sapply(c("A", "B", "C"), function(x) rep(x, 10))),
  as.data.frame(matrix(rnorm(3000, mean=7, sd=0.25), nrow=30))
)
colnames(d) <- c("sample", "sample_group", paste0("gene", 1:100))
rownames(d) <- d$sample

# Complete the code in the following function to run ANOVA on the gene expression data generated above.

runAnova <- function(input_data) {
  Gene<-data.frame(colnames(input_data[3:ncol(input_data)]))
  GroupA <- input_data[input_data$sample_group == "A",]
  print(str(Gene))
  print(str(GroupA))
  n<-c()
  P<-c()
  P_adj<-c()
  for (i in (3:length(colnames(input_data)))){
    samp<-length(GroupA[,i])
    n<-c(n,samp)
    res_aov <- aov(input_data[,i] ~ input_data$sample_group, data = input_data )
    p_val <- summary(res_aov)[[1]][[1,"Pr(>F)"]]
    adj<-p.adjust(p_val, method = "BH", n=samp)
    P <-c(P,p_val)
    P_adj<-c(P_adj,adj)
  }
  data<-cbind(Gene,n,P,P_adj)
  colnames(data)<-c("Gene","n","P","P_adj")
  return(data)
  # TODO: write the codes to compute the ANOVA test on the gene expression data generated in the format above. Run the ANOVA test for each of the gene against the sample group column.The function should return a data.frame of the results in with the column names "Gene", "n", "P" and "P_adj" which represent the gene name, number of samples, P value and Benjamini and Hochberg corrected P values. (Q12 6 marks)
}
write.table(runAnova(d), file="q12.tsv", quote=F, na="", sep="\t", col.names=T, row.names=F)

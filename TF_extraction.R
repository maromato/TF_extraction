#Initially, you need to normalize raw microarray data and make a spread sheet for gene expression as shown elsewhere.

#input spread sheet for microarray data

x = read.table("XXXX.txt",header=T,sep="\t")

x2 = x[,2:ncol(x)]

x2 = as.matrix (x2)

rownames(x2)

#statistical analysis for differential expression

library(limma)
design = cbind(Cell_A = c(1,1,1,0,0,0),Cell_B = c(0,0,0,1,1,1))

fit = lmFit(x2, design=design)# Fit the original matrix to the above design.

contrastsMatrix = makeContrasts("Cell_B-Cell_A",
levels = design)# We want to compare A vs. B, A vs. C and B vs. C

fit2 = contrasts.fit(fit, contrasts = contrastsMatrix) # Making the comparisons.

out = eBayes(fit2) # Moderating the t-tetst by eBayes method.

p.value = out$p.value #to put p.values for indevisual genes into the vector, p.value

q.value = apply(p.value, MARGIN=2, p.adjust, method="BH")#to put
q.values for indevisual genes into the vector, a.value

ranking = apply(p.value, MARGIN=2, rank)#to put the ranking for indevisual genes in term of p.value into the vector, ranking
tmp = cbind(x2, p.value, q.value, ranking)

#extract mouse TF data from original spread sheet.

b = read.table("mouseTF.txt",sep=" \t",header=F)#you need to make a entire list of the mouse transcription factors from the information
in the site (http://genome.gsc.riken.jp/TFdb/tf_list.html)

obj = rownames(x2) %in% b[,1] # to find the rows that are included in the list of transcription factors (https://www.biostars.org/p/5737/)

tmp2 = tmp[obj,]

ob = ranking < 65

tmp3 = tmp2[ob,]

z = tmp3[,1:6]

#Clustering and visualization.

library(gplots)

heatmap.2(as.matrix(z), col=greenred(75), scale="row", key=T, keysize=1.5,density.info="none", trace="none",cexCol=0.9, cexRow=0.5)

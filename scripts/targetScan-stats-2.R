#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly=T)
olddir <- getwd()
data.path <- args[1]
n <- as.integer(args[2])
percent <- as.integer(args[3])
setwd(data.path)
print(paste("Getting in in ", getwd(), sep=""))
lf <- list.files(pattern=".data$")

#Ploting %L genes & %S genes and -log10(p-value) of P(X >= x; N, K, n)

probs_l <- c()
probs_s <- c()
i <- c(0)
selected.mirnas <- c()
matrix_percent <- matrix(rep(0, 2*length(lf)), nrow=2, ncol=length(lf))
print(length(lf))
str(matrix_percent)
for (file in lf) {
	table <- read.table(file, header=F)
	print("Table read")
	sample <- table[,1]
	K <- c(4578) # complete dataset
#	K <- length(grep(".S$", sample))  ### Change if L or S
#	K <- c(12951) # all dataset
#	K <- c(2967) # singletons dataset
#	K <- c(5614) # ONLY ONLY ONLY PARALOGS
	str(K)
	N <- c(9156) # complete dataset
#	N <- length(sample) 
#	N <- c(22530) # all dataset
#	N <- c(3761) # singletons dataset
#	N <- c(11228) # ONLY ONLY ONLY PARALOGS
	str(N)
	subsample_size <- (percent / 100) * length(sample)
	subsample <- sample[1:subsample_size]
	x_l <- length(grep(".L$", subsample)) ### Change if L or S

# Change if LS
	x_s <- length(grep(".S$", subsample))
#	x_x <- length(subsample) - x_l
	prb_l <- phyper(x_l, K, N-K, subsample_size, lower.tail=F)

# Change if LS
	prb_s <- phyper(x_s, N-K, K, subsample_size, lower.tail=F)
#	prb_s <- phyper(x_x, N-K, K, subsample_size, lower.tail=F)
	i = i + 1
	if (-log10(prb_s) * 5 < -log10(prb_l)) {
		selected.mirnas <- c(selected.mirnas, file)
	}
	matrix_percent[1,i] <- (x_l/subsample_size) * 100

# Chenge if LS
	matrix_percent[2,i] <- (x_s/subsample_size) * 100
#	matrix_percent[2,i] <- (x_x/subsample_size) * 100
	probs_l <- c(probs_l, -log10(prb_l))
	probs_s <- c(probs_s, -log10(prb_s))
	x_axis <- 1:(length(lf))
}
# probs and matrix_percent are filled and i is equal to N

# Change if LS
rownames(matrix_percent) <- c("%_L_genes", "%_S_genes")
#rownames(matrix_percent) <- c("%_L_genes", "%_X_genes")
colnames(matrix_percent) <- lf
print(matrix_percent)
print(selected.mirnas)

pdf("stackedBars.pdf", width=10, height=5)
par(mar=c(3,5,3,5))
barplot(as.matrix(matrix_percent), axes=F, ylim = c(0, 100), xlab="", ylab="")
axis(2, ylim=c(0, 100), col = "black", las=1)
mtext("% of L & S genes", side=2, line=2.5)
box()
par(new=T)
plot(x_axis, probs_l, col = "blue", ylim = c(0, max(  max(probs_l[!is.infinite(probs_l)]) , max(probs_s[!is.infinite(probs_s)]) )), ylab="", xlab="", type="l", axes=F)
lines(x_axis, probs_s, col="red", ylim = c(0, max(  max(probs_l[!is.infinite(probs_l)]) , max(probs_s[!is.infinite(probs_s)]) )), ylab="", xlab="", axes=F)
axis(4, ylim=c(0, max( max(probs_l[!is.infinite(probs_l)]) , max(probs_s[!is.infinite(probs_s)]) )), ylab="", xlab="", col = "blue", las=1)
abline(h=1.4, col="red")
mtext("-log10(p-value); P(X >= x; N,K,n)", side=4, line=2.5, col="blue")
#axis(1, pretty(1:length(lf)))
mtext("miRNAs family", side=1, line=2.5)
dev.off()


######################################################################

#//
#	Ploting INtersection of target paralogs and significance of the intersection.
#//


######################################################################

columns <- c(5)
total <- length(lf)
rows <- total/columns
pdf("targetScan.pdf", width=columns*4, height=rows*3)
par(mfrow=c(rows, columns))
par(mar=c(3,5,3,5))
for (file in lf) {
	print(paste("Accessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene", "cs++")
	rownames(table) <- 1:length(table[,1])
	array_of_genes <- as.character(table$gene)
	N <- length(array_of_genes)
	print("N")
	print(N)
	window <- n #<- N/200
	print("k")
	k <- length(grep(".L", array_of_genes))
	#k <- N/2
	print(k)
	x_axis <- c()
	y_l_axis <- c()
	y_s_axis <- c()
	y2_axis <- c()
	#j <- c(1)
	for (i in seq(1, N-n, by=n)) {
		subarray_of_genes <- array_of_genes[1:window]
		print(head(subarray_of_genes))
		print("Len subarray")
		print(length(subarray_of_genes))
		x_l <- length(grep(".L$", subarray_of_genes))
		print("Total L in sample")
		print(x_l)
		x_s <- length(grep(".S$", subarray_of_genes))
		print("prb_l")
		prb_l <- phyper(x_l, k, N-k, window, lower.tail=F)
		print(prb_l)
		#print(paste("y_l_axis[j] ", y_l_axis[j], " + prb_l ", prb_l, " = ", y_l_axis[j] + prb_l ))
		#prb_l <- y_l_axis[j] + prb_l
		prb_s <- phyper(x_s, N-k, k, window, lower.tail=F)
		#prb_s <- y_s_axis[j] + prb_s
		x_axis <- c(x_axis, window -  (n/2))
		y_l_axis <- c(y_l_axis, -log10(prb_l))
		y_s_axis <- c(y_s_axis, -log10(prb_s))
		y2_axis <- c(y2_axis, (x_l/window)*100)
		window = window + n
		#j = j + 1
	}
	min_y <- min(min(y_l_axis[!is.infinite(y_l_axis)]), min(y_s_axis[!is.infinite(y_s_axis)]))
	print(min_y)
	max_y <- max(max(y_l_axis[!is.infinite(y_l_axis)]), max(y_s_axis[!is.infinite(y_s_axis)]))
	print(max_y)
	min_y2 <- min(y2_axis)
	max_y2 <- max(y2_axis)
	#print(x_axis)
	#print(y_axis)
	q1 <- (25*N)/100
	mode <- (50*N)/100
	q2 <- (75*N)/100
	plot(x_axis, y_l_axis, main=file, ylim=c(min_y, max_y), xlab="", ylab="", type="l", axes=F)
	lines(x_axis, y_s_axis, ylim=c(min_y, max_y), xlab="", ylab="", axes=F, col="gray")
	axis(2, ylim=c(min_y, max_y), col="black", las=1)
	mtext("P(X <= x)", side=2, line=2.5)
	box()
	par(new=T)
	plot(x_axis, y2_axis, col="blue", ylim=c(min_y2,max_y2), xlab="", ylab="", axes=F)
	axis(4, ylim=c(min_y2, max_y2), col="blue", las=1)
	mtext("Percent of L genes", side=4, line=2.5)
	axis(1, pretty(1:N))
	mtext("Genes", side=1, line=2.5, col="black")
	abline(v=q1, col="red")
	abline(v=mode, col="red")
	abline(v=q2, col="red")
}
dev.off()

.silenced <- function (){
pdf("targetScan_boxplot.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for (file in lf) {
	print(paste("Accessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene", "cs++")
	rownames(table) <- 1:length(table[,1])
	array_of_genes <- as.character(table$gene)
	for (i in c(0, 25, 50, 75)) {
		subarray_of_genes <- array_of_genes[i:i+25]
		L <- length(grep(".L", subarray_of_genes))
		S <- length(grep(".S", subarray_of_genes))

	}
	#print(x_axis)
	#print(y_axis)
	q1 <- (25*N)/100
	mode <- (50*N)/100
	q2 <- (75*N)/100
	plot(x_axis, y_axis, main=file, xlab="Genes", ylab="P(X <= x)", type="l")
	abline(v=q1, col="red")
	abline(v=mode, col="red")
	abline(v=q2, col="red")
}
dev.off()
} #END .silenced


lf <- list.files(pattern=".sorted$")
pdf("targetScan_paralogs.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for (file in lf) {
	print(paste("Accsessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene_name", "geneL", "csL", "geneS", "csS")
	rownames(table) <- 1:length(table[,1])
	table <- table[order(table$csL),]
	x_axis <- 1:length(table[,1])
	y_axis <- table$csS
	plot(x_axis, y_axis, main=file, xlab="L genes", ylab="CS++ for S genes")
	abline(h=median(y_axis), col="red")
	abline(h=median(table$csL), col="blue")
}
dev.off()

lf <- list.files(pattern=".sorted$")
pdf("targetScan_paralogs_boxplot.pdf", width=columns*3, height=rows*3)
par(mfrow=c(rows, columns))
for (file in lf) {
	print(paste("Accsessing ", file, sep=""))
	table <- read.table(file, header=F)
	colnames(table) <- c("gene_name", "geneL", "csL", "geneS", "csS")
	rownames(table) <- 1:length(table[,1])
	table <- table[order(table$csL),]
	N <- length(table[,1])
	q1 <- (25*N)/100
	median <- (50*N)/100
	q2 <- (75*N)/100
	q1.table <- table[0:q1,]
	median.table <- table[q1:median,]
	q2.table <- table[median:q2,]
	last.table <- table[q2:N,]
	y_min <- min(min(table$csL), min(table$csS))
	y_max <- max(max(table$csL), max(table$csS))
	boxplot(q1.table$csL, main=file, at=1, xlim=c(0,8), ylim=c(y_min - 1, y_max + 1), col="red")
	boxplot(q2.table$csS, at=2, add=T, col="red")
	boxplot(median.table$csL, at=3, add=T, col="blue")
	boxplot(median.table$csS, at=4, add=T, col="blue")
	boxplot(q2.table$csL, at=5, add=T, col="violet")
	boxplot(q2.table$csS, at=6, add=T, col="violet")
	boxplot(last.table$csL, at=7, add=T, col="green")
	boxplot(last.table$csS, at=8, add=T, col="green")
}
dev.off()

setwd(olddir)
print(paste("Returning to ", getwd(), sep=""))



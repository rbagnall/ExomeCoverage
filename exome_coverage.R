#!/bin/env Rscript

# plot histogram of per base coverage across the exome, and sex assignment

args <- commandArgs(TRUE)
verbose = TRUE

coverageFile = args[1]
# read per.base.coverage as a table (/path/to/sample.per.base.coverage)
sample_name <- unlist(strsplit(args[1], split=".per.base.coverage.gz", fixed=TRUE)) # name of /path/to/sample

sample_id <- unlist(strsplit(sample_name, split="/", fixed=TRUE))
data <- read.table(gzfile(args[1]), colClasses=c("character", rep("integer", 4)))
# read in coverage file
coverage <- data[,5]
# get maximum coverage value
max_num <-max(coverage)
# get number of bases in the exome
covlength <- length(coverage)
# get number of bases with 0 coverage
zerocoverage = length(which(coverage==0))
# get number of bases with less than 10 reads
lessthantenbases = length(which(coverage<10))
# get number of bases with less than 10 reads
lessthantwentybases = length(which(coverage<20))
# get number of bases with less than 50 reads
lessthanfiftybases = length(which(coverage<50))
# get maximum frequency value
histinfo <-hist(coverage, plot=FALSE, breaks=max_num)
# maximum frequency
max_frequency <-max(histinfo$counts)
# covearge on X chromosome
chrX <- data[grep("^X", data[,1]), ]
# average coverage on X chromosome
aveX <- round(mean(chrX[,5]))
# covearge on Y chromosome
chrY <- data[grep("^Y", data[,1]), ]
# average coverage on Y chromosome
aveY <- round(mean(chrY[,5]))
# gender ratio
sexratio <- paste("(", aveX, ":", aveY, ")", sep="")
# calculate gender - if ratio of coverage on Y is less than 0.3 of X, then its a female
if (aveY/aveX < 0.3) {
    gender <- "Female    "
} else {
    gender <- "Male      "
}
# calculate X axis position to print stats
plotXposition = round(mean(coverage*3))*0.9

# calculate Y axis positions to print stats
plotYgender = max_frequency*0.65
plotYone = max_frequency*0.6
plotYtwo = max_frequency*0.55
plotYthree = max_frequency*0.5
plotYfour = max_frequency*0.45
plotYfive = max_frequency*0.4
plotYsix = max_frequency*0.35
# calculate % bases having 0 reads
zerocov <- round((zerocoverage/covlength)*100, digits=1)
# calculate % bases having 10 or more reads
readabletencoverage = round((1-(lessthantenbases/covlength))*100, digits=1)
# calculate % bases having 20 or more reads
readabletwentycoverage = round((1-(lessthantwentybases/covlength))*100, digits=1)
# calculate % bases having 50 or more reads
readablefiftycoverage = round((1-(lessthanfiftybases/covlength))*100, digits=1)

# create .pdf file to print to
outfile= paste(coverageFile, ".pdf", sep="")
# create .txt file to print to
outtxt <- file((paste(sample_name, ".coverage.metrics.txt", sep="")))
writeLines(c(paste("Sample Name", "Sex", "Y:X coverage", "Mean Coverage", "Median Coverage", "Max Coverage", "Depth=0 (%)", "Depth >= 10 reads (%)", "Depth >= 20 reads (%)", "Depth >=50 read (%)", sep="\t"), paste(tail(sample_id, n=1), gender, format(aveY/aveX, digits=2, nsmall=2), round(mean(coverage),1), round(median(coverage),1), max(coverage), zerocov, readabletencoverage, readabletwentycoverage, readablefiftycoverage, sep="\t")), outtxt)
close(outtxt)
pdf(outfile)
# plot histogram
par(mar=c(5,6,4,2))
hist(coverage, col="lightblue", ylim=c(0,(max_frequency*1.2)),xlim=c(0,(round(mean(coverage*3)))), breaks=max_num, las=1, cex.axis=1, border="lightblue", xlab="Per base sequence read depth", ylab="", main=paste(tail(sample_id, n=1), " Histogram of Coverage", sep=""))
# add Y axis label
mtext("Number of bases at read depth", 2, line=5)
# print stats onto the histogram
text(plotXposition, plotYgender, gender, cex=0.9, pos=2)
text(plotXposition, plotYgender, sexratio, cex=0.9)
text(plotXposition, plotYone, "Mean coverage =  ", cex=0.9, pos=2)
text(plotXposition, plotYone, round(mean(coverage),1), cex=0.9)
text(plotXposition, plotYtwo, "Median coverage =  ", cex=0.9, pos=2)
text(plotXposition, plotYtwo, round(median(coverage),1), cex=0.9)
text(plotXposition, plotYthree, "Max coverage =  ", cex=0.9, pos=2)
text(plotXposition, plotYthree, max(coverage), cex=0.9)
text(plotXposition, plotYfour, "Depth >=10 reads (%) =  ", cex=0.9, pos=2)
text(plotXposition, plotYfour, readabletencoverage, cex=0.9)
text(plotXposition, plotYfive, "Depth >=20 reads (%) =  ", cex=0.9, pos=2)
text(plotXposition, plotYfive, readabletwentycoverage, cex=0.9)
text(plotXposition, plotYsix, "Depth >=50 reads (%) =  ", cex=0.9, pos=2)
text(plotXposition, plotYsix, readablefiftycoverage, cex=0.9)

dev.off()
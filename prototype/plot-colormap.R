# Clear workspace
rm(list=ls())

# Strings as factors
options(stringsAsFactors = F)

# Import
library(optparse, quietly=TRUE)

# Read arguments
option_list <- list(
make_option(c("-d","--data"), help = "Color map (*.colormap)", metavar="FILE"),
make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."))

parsed <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

if(!exists("data", where=parsed)) {
    print("Error: Invalid arguments. Check help (-h).")
    quit("no", 1)
}

# Read data
if(parsed$verbose) print(paste("Reading file", parsed$data))

cm_data <- read.csv(parsed$data, header=TRUE, comment.char='#', na.strings="NA")
cm_data <- cm_data[complete.cases(cm_data),]
outf <- paste(gsub("\\.", "-", parsed$data), ".pdf", sep="")

# Plot colormap
if(parsed$verbose) print("Plotting colormap")

pdf(outf)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, legend=as.character(cm_data$value), fill=as.character(cm_data$color), xjust=0.5, yjust=0.5)
title(parsed$data)

# Write out
junk <- dev.off()
if(parsed$verbose) print(paste("Wrote file:", outf))

# Warn
wa <- warnings()
if(class(wa) != "NULL")
    print(wa)


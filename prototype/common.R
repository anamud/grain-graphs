# Libraries
library(ggplot2, quietly=TRUE) # Plotting
library(reshape2, quietly=TRUE) # Melting
library(scales, quietly=TRUE) # Log transformations
library(RColorBrewer, quietly=TRUE) # Colors
suppressMessages(library(igraph, quietly=TRUE))
suppressMessages(library(gdata, quietly=TRUE, warn.conflicts=FALSE))
library(optparse, quietly=TRUE)

# Treat warnings as errors
#options(warn=2)

# Set number of precision digits
options("scipen"=999)

# Strings as factors
options(stringsAsFactors = F)

# Timing functions
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
    type <- match.arg(type)
    assign(".type", type, envir=baseenv())
    if(gcFirst) gc(FALSE)
    tic <- proc.time()[type]
    assign(".tic", tic, envir=baseenv())
    invisible(tic)
}

toc <- function(message)
{
    type <- get(".type", envir=baseenv())
    toc <- proc.time()[type]
    tic <- get(".tic", envir=baseenv())
    cat(sprintf("%s: %f sec\n", message, toc - tic))
    invisible(toc)
}

# Useful functions
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
my_print <- function(x="") {cat(paste(x, "\n", sep=""))}

# Plotting helpers
common_theme <- theme(
                      text = element_text(size=18),
                      axis.text.x = element_text(colour = "black",
                                                 size=18,
                                                 angle = 0,
                                                 hjust = 0),
                      axis.text.y = element_text(colour = "black",
                                                 size=18),
                      legend.position="top",
                      legend.title = element_blank(),
                      legend.text = element_text(size=18),
                      plot.background = element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour="black",
                                                  fill=NA)
                      )

bar_plotter <- function(df, log=F, xt, yt, mt=" ", tilt=F, tilt_angle=0)
{
    colnames(df) <- c("subject", "variable")
    df.melt <- melt(df, id.vars="subject")

    if(tilt == T)
    {
        barp_t <- common_theme + theme(axis.text.x = element_text(angle = tilt_angle, hjust = 1))
    } else {
        barp_t <- common_theme
    }

    if(length(df$subject) > 20) barp_t <- barp_t + theme(axis.text.x = element_text(size=18))

    barp <- ggplot(df.melt[complete.cases(df.melt),], aes(x=subject,y=value)) +
        geom_bar(stat="identity", width = 0.8, position = position_dodge(width = 0.8))
    if(log == T) {
        barp_l <- labs(x=xt, y=paste(yt,", log scale", sep=""),title=mt)
        barp <- barp + scale_y_log10(labels=scientific) + barp_l + barp_t
    } else {
        barp_l <- labs(x=xt, y=yt,title=mt)
        barp <- barp + scale_y_continuous(labels=scientific) + barp_l + barp_t
    }

    print(barp)
}

bar_plotter_two_bars <- function(df, xt, yt, mt=" ", lt=c("variable1", "variable2"), tilt=F, tilt_angle=0, log=F)
{
    colnames(df) <- c("subject", lt)
    df.melt <- melt(df, id.vars="subject")

    if(tilt == T) {
        barp_t <- common_theme + theme(axis.text.x = element_text(angle = tilt_angle, hjust = 1))
    } else {
        barp_t <- common_theme
    }

    if(length(df$subject) > 20) barp_t <- barp_t + theme(axis.text.x = element_text(size=8))

    barp <- ggplot(df.melt[complete.cases(df.melt),], aes(x=subject, y=value, fill=factor(variable))) +
        geom_bar(stat="identity", width = 0.8, position = position_dodge(width = 0.8)) +
        #scale_fill_brewer(palette="Dark2")
        scale_fill_manual(values=c("black","gray80"))

    if(log == T) {
        barp_l <- labs(x=xt, y=paste(yt,", log scale", sep=""),title=mt)
        barp <- barp + scale_y_log10(labels=scientific) + barp_l + barp_t
    } else {
        barp_l <- labs(x=xt, y=yt,title=mt)
        barp <- barp + scale_y_continuous(labels=scientific) + barp_l + barp_t
    }

    print(barp)
}

box_plotter <- function(v, xt, yt, mt=" ",log=F)
{
    v_comp <- v[complete.cases(v)]
    v_comp <- v_comp[is.finite(v_comp)]
    v_comp <- v_comp[!is.nan(v_comp)]

    if(length(v_comp) > 0) {
        v.melt <- melt(v)
        if(log == T) {
            boxp <- ggplot(v.melt, aes(y=value,x=factor(0))) +
                #geom_boxplot(outlier.size=0) +
                #geom_jitter(position=position_jitter(width=.2), size=2) +
                geom_boxplot() +
                labs(y=paste(yt,"log scale",sep=", "),x=xt,title=mt) +
                scale_y_log10(labels=scientific) +
                common_theme + theme(axis.text.x = element_blank())
        } else {
            boxp <- ggplot(v.melt, aes(y=value,x=factor(0))) +
                #geom_boxplot(outlier.size=0) +
                #geom_jitter(position=position_jitter(width=.2), size=2) +
                geom_boxplot() +
                labs(y=yt,x=xt,title=mt) +
                scale_y_continuous(labels=scientific) +
                common_theme + theme(axis.text.x = element_blank())
        }

        print(boxp)
    } else {
        print("Cannot plot. No complete cases.")
    }
}

box_plotter_two_boxes <- function(v1, v2, xt, yt, mt=" ", lt=c("v1", "v2"), log=F)
{
    v1_comp <- v1[complete.cases(v1)]
    v1_comp <- v1_comp[is.finite(v1_comp)]
    v1_comp <- v1_comp[!is.nan(v1_comp)]

    v2_comp <- v2[complete.cases(v2)]
    v2_comp <- v2_comp[is.finite(v2_comp)]
    v2_comp <- v2_comp[!is.nan(v2_comp)]

    if(length(v1_comp) > 0 & length(v2_comp) > 0) {
        v1_comp <- data.frame(group=lt[1], value=v1_comp)
        v2_comp <- data.frame(group=lt[2], value=v2_comp)
        v12 <- rbind(v1_comp, v2_comp)
        v12.melt <- v12

        if(log == T) {
            boxp <- ggplot(v12.melt, aes(y=value,x=factor(group))) +
                geom_boxplot() +
                labs(y=paste(yt,"log scale",sep=", "),x=xt,title=mt) +
                scale_y_log10(labels=scientific) +
                common_theme
        } else {
            boxp <- ggplot(v12.melt, aes(y=value,x=factor(group))) +
                geom_boxplot() +
                labs(y=yt,x=xt,title=mt) +
                scale_y_continuous(labels=scientific) +
                common_theme
        }

        print(boxp)
    } else {
        print("Cannot plot. No complete cases.")
    }
}

# Modify alpha color attribute (A) in RGBA
# From http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add.alpha <- function(col, alpha=1){
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
}


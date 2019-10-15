# Headers
source("src/headers.R")

# Wrapper functions
.combined.plot <- function(quantile, dims=1:2, clean=5, abs=FALSE){
    .pca.plot <- function(pca, dims, scale=5, ...){
        density <- kde(pca$x[,dims])
        plot(pca$x[,dims], type="n", axes=FALSE)
        axis(1)
        axis(2)
        plot(density, display="filled.contour", cont=seq(10,90,10), cols=heat.colors(length(seq(10,90,10))), add=TRUE, axes=FALSE)
        points(pca$x, pch=20, cex=.5, col="grey60")
        arrows(rep(0,ncol(pca$x)),rep(0,ncol(pca$x)), pca$rotation[,dims[1]]*scale, pca$rotation[,dims[2]]*scale, length=.1, lwd=3)
        #arrows(rep(0,ncol(pca$x)),rep(0,ncol(pca$x)), -pca$rotation[,dims[1]]*scale, -pca$rotation[,dims[2]]*scale, length=.1, lwd=3)
        text.locs <- pca$rotation*scale*1.3
        text(text.locs[,dims[1]], text.locs[,dims[2]], names(pca$center), font=2)
    }

    index <- .simplify.group("bootstrap.index", quantile, clean=clean,  abs=abs)[,1:9]
    present <- .simplify.group("present", quantile, abs=abs)[,1:9]
    past <- .simplify.group("past", quantile, abs=abs)[,1:9]
    
    layout(matrix(1:6, nrow=2, byrow=TRUE), widths=c(5,5,5), heights=c(5,1))
    par(mar=c(1,1,1,1))
    present <- prcomp(present, scale=TRUE)
    past <- prcomp(past, scale=TRUE)
    index <- prcomp(index, scale=TRUE)
    .pca.plot(past, 1:2, bty="o")
    .pca.plot(present, 1:2)
    .pca.plot(index, 1:2)
    plot(past, main="")
    plot(past, main="")
    plot(index, main="")
}

# Load data
plants <- readRDS("clean-data/plants-index.RDS")
fungi <- readRDS("clean-data/plants-index.RDS")
insects <- readRDS("clean-data/plants-index.RDS")
mammals <- readRDS("clean-data/mammals-index.RDS")
reptiles <- readRDS("clean-data/reptiles-index.RDS")
birds <- readRDS("clean-data/birds-index.RDS")

# PCAs
pdf("figures/pca-05.pdf", width=15, height=6); .combined.plot("0.05"); dev.off()
pdf("figures/pca-25.pdf", width=15, height=6); .combined.plot("0.25"); dev.off()
pdf("figures/pca-50.pdf", width=15, height=6); .combined.plot("0.5"); dev.off()
pdf("figures/pca-75.pdf", width=15, height=6); .combined.plot("0.75"); dev.off()
pdf("figures/pca-95.pdf", width=15, height=6); .combined.plot("0.95"); dev.off()

# Simulating PCAs
correls <- numeric(100)
real <- prcomp(combined.pres[,1:9], scale=TRUE)
for(i in seq_along(correls)){
    past <- combined.past[,1:9]
    for(i in seq_len(9))
        past[,i] <- past[,i] + rnorm(nrow(past), sd=sd(abs(combined.pres[,i] - combined.past[,i])))
    pres <- past
    for(i in seq_len(9))
        past[,i] <- past[,i] + rnorm(nrow(past), sd=sd(abs(combined.pres[,i] - combined.past[,i])))
    diff <- pres - past
    past.pca <- prcomp(past, scale=TRUE)
    diff.pca <- prcomp(diff, scale=TRUE)
    correls[i] <- mantel.test(abs(diff.pca$rotation), abs(past.pca$rotation))$z.stat
}

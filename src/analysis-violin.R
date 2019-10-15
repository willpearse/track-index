# Headers
source("src/headers.R")

# Load data
plants <- readRDS("clean-data/plants-index.RDS")
fungi <- readRDS("clean-data/plants-index.RDS")
insects <- readRDS("clean-data/plants-index.RDS")
mammals <- readRDS("clean-data/mammals-index.RDS")
reptiles <- readRDS("clean-data/reptiles-index.RDS")
birds <- readRDS("clean-data/birds-index.RDS")

# Plots
pdf("figures/violin-tmp-50.pdf")
index <- .simplify.group("index", quantile="0.5", clean=10, abs=FALSE)
se <- .simplify.group("mad.index", quantile="0.5", abs=FALSE)
se <- se[se$tmp < 5,]
se <- se[rownames(se) %in% rownames(index),]
index <- index[rownames(index) %in% rownames(se),]
data <- data.frame(index=index$tmp, se=se$tmp, taxon=index$taxon)
cols <- setNames(c("forestgreen", "black", "blue", "orange", "red", "grey60"), unique(index$taxon))
data$col <- alpha(cols[data$taxon], alpha=1-((data$se^2)/max(data$se^2)))
data <- data[order(data$taxon, data$index),]
data$x <- as.numeric(factor(data$taxon))
counts <- table(data$taxon)
curr <- 0
for(i in seq_along(counts)){
    taxon <- names(counts)[i]
    data$x[data$taxon==taxon] <- data$x[data$taxon==taxon] + seq(-.3, .3, length.out=counts[i])
}
medians <- with(data, tapply(index, taxon, quantile, probs=.5))
lower <- with(data, tapply(index, taxon, quantile, probs=.25))
upper <- with(data, tapply(index, taxon, quantile, probs=.75))
with(data, plot(index ~ x, xlab="", ylab=expression(track[50]~temperature), axes=FALSE, pch=20, col=col))
arrows(1:6-.2, medians, 1:6+.2, col=cols, length=0, lwd=5)
arrows(1:6-.2, lower, 1:6+.2, col=cols, length=0, lwd=5)
arrows(1:6-.2, upper, 1:6+.2, col=cols, length=0, lwd=5)
abline(h=0, col="grey40", lwd=3, lty=2)
abline(h=1, col="grey40", lwd=3, lty=2)
axis(1, labels=names(cols), at=1:6)    
axis(2)
dev.off()

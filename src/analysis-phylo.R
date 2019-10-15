# Headers
source("src/headers.R")

# Wrapper functions
.make.comparative.data <- function(data, tree, merge=TRUE){
    tree$tip.label <- tolower(gsub(" ", "_", tree$tip.label))
    tree$node.label <- NULL
    data <- as.data.frame(data)
    data$species <- rownames(data)
    data$taxon <- NULL
    data <- data[!duplicated(data$species),]
    if(merge)
        tree <- congeneric.merge(tree, data$species)
    return(comparative.data(tree, data, species))
}
.phy.sig <- function(c.data, method="lambda")
    sapply(c.data$data, function(x) phylosig(c.data$phy, x, method))
.dom.niche <- function(c.data)
    setNames(
        apply(as.matrix(c.data$data), 1, function(x) which.min(abs(x))),
        rownames(c.data$data)
    )
.phy.sig.wrap <- function(data, tree, method="lambda", merge=TRUE, clean=5){
    output <- matrix(NA, ncol=9, nrow=5)
    rownames(output) <- dimnames(data)[[2]]
    colnames(output) <- dimnames(data)[[1]]
    for(i in seq_len(nrow(output))){
        c.data <- .make.comparative.data(.simplify(data, "bootstrap.index", rownames(output)[i], clean=clean), tree, merge)
        output[i,] <- .phy.sig(c.data)[1,]
    }
    return(output)
}
.dom.niche.wrap <- function(data, tree, method="lambda", merge=TRUE, clean=5){
    output <- setNames(numeric(nrow(data)), dimnames(data)[[1]])
    for(i in seq_len(nrow(output))){
        c.data <- .make.comparative.data(.simplify(data, "bootstrap.index", rownames(output)[i], clean=clean), tree, merge)
        model <- fitDiscrete(multi2di(c.data$phy), model="SYM", .dom.niche(c.data), transform="lambda")
        output[i] <- model$opt$lambda
    }
    return(output)
}
.get.branch.trait <- function(trait, tree){
    trait <- c(trait[tree$tip.label], fastAnc(tree, trait))
    names(trait)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
    edge.mat <- matrix(trait[tree$edge], nrow(tree$edge), 2)
    trait <- rowMeans(edge.mat)
    return(trait)
}
.limit <- function(x){
    min <- min(x); max <- max(x)
    x <- x - min
    x <- x/max(x)
    attr(x, "range") <- c(min, max)
    return(x)
}

# Load data
plants <- readRDS("clean-data/plants-index.RDS")
fungi <- readRDS("clean-data/plants-index.RDS")
insects <- readRDS("clean-data/plants-index.RDS")
mammals <- readRDS("clean-data/mammals-index.RDS")
reptiles <- readRDS("clean-data/reptiles-index.RDS")
birds <- readRDS("clean-data/birds-index.RDS")

# Phylogeny
c.plant.index <- .make.comparative.data(.simplify(plants, "bootstrap.index", "0.5", clean=5), read.tree("phylogenies/ALLOTB.tre"), FALSE)
c.plant.past <- .make.comparative.data(.simplify(plants, "past", "0.5"), c.plant.index$phy, FALSE)



past.edge <- .get.branch.trait(setNames(c.plant.past$data$tmp,rownames(c.plant.past$data)), c.plant.past$phy)
cut.past.edge <- as.numeric(cut(.limit(past.edge), breaks=seq(-.001,1.001,.001)))
cols <- hcl.colors(length(seq(-.001,1.001,.001)))

setwd("figures")
pdf("phylo-tmp-50.pdf")
plot(c.plant.index$phy, type="fan", show.tip.label=FALSE, no.margin=TRUE, edge.color=cols[cut.past.edge])
cut.index <- as.numeric(cut(.limit(c.plant.index$data$tmp), breaks=seq(-.001,1.001,.001)))
willeerd.tiplabels(pch=20, col=cols[cut.index], radial.adj=1.02)
dev.off()

print("plants")
print(.phy.sig.wrap(plants, read.tree("phylogenies/ALLOTB.tre")))
print("mammals")
print(.phy.sig.wrap(plants, read.tree("phylogenies/faurby_mammals.tre"))print(
print("reptiles")
print(.phy.sig.wrap(plants, read.tree("phylogenies/squamate_tree.tre")))


# Find the dominant trait
library(geiger)
fitDiscrete(multi2di(c.mammal$phy), .dom.niche(c.mammal), transform="lambda")
fitDiscrete(multi2di(c.mammal$phy), model="SYM", .dom.niche(c.mammal), transform="lambda")

phylosig(c.plant$phy, .dom.niche(c.plant), method="K")

# Find dominant per species
.dom.niche(c.plant)

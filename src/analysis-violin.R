# Headers
source("src/headers.R")

# Load data

plants <- readRDS("clean-data/plants-index.RDS")
fungi <- readRDS("clean-data/plants-index.RDS")
insects <- readRDS("clean-data/plants-index.RDS")
mammals <- readRDS("clean-data/mammals-index.RDS")
reptiles <- readRDS("clean-data/reptiles-index.RDS")
birds <- readRDS("clean-data/birds-index.RDS")

# Do the plots
.violin.plot <- function(climate, quant, clean=10, se.thresh=5, abs=FALSE, index="bootstrap.index"){
    # Load and match data
    index <- .simplify.group(index, quantile=quant, clean=clean, abs=abs)
    se <- .simplify.group("mad.index", quantile=quant, abs=abs)
    se <- se[se[,climate] < se.thresh,]
    se <- se[rownames(se) %in% rownames(index),]
    index <- index[rownames(index) %in% rownames(se),]
    data <- data.frame(index=index[,climate], se=se[,climate], taxon=index$taxon)

    # Prepare colours and order for plotting
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

    # Lookups for nicer names
    c.var <- setNames(c("cloud-cover","frost-days","potential-evapotranspiration","precipitation","min-temperature","mean-temperature","max-temperature","vapour-pressure"),
                      c("cld","frs","pet","pre","tmn","tmp","tmx","vap")
                      )[climate]
    q.var <- setNames(c("5th","25th","50th","75th","95th"), c("0.05","0.25","0.5","0.75","0.95"))[quant]

    # Calculate plot-points
    medians <- with(data, tapply(index, taxon, quantile, probs=.5))
    lower <- with(data, tapply(index, taxon, quantile, probs=.25))
    upper <- with(data, tapply(index, taxon, quantile, probs=.75))

    # Do work
    eval(substitute(
        with(data, plot(index ~ x, xlab="", ylab=expression(track[XXX](YYY)), axes=FALSE, pch=20, col=col)),
        list(XXX=q.var, YYY=c.var)
    ))
    arrows(1:6-.2, medians, 1:6+.2, col=cols, length=0, lwd=5)
    arrows(1:6-.2, lower, 1:6+.2, col=cols, length=0, lwd=5)
    arrows(1:6-.2, upper, 1:6+.2, col=cols, length=0, lwd=5)
    abline(h=0, col="grey40", lwd=3, lty=2)
    abline(h=1, col="grey40", lwd=3, lty=2)
    axis(1, labels=names(cols), at=1:6)    
    axis(2)
}
for(climate in c("cld","frs","pet","pre","tmn","tmp","tmx","vap")){
    for(quant in c("0.05","0.25","0.5","0.75","0.95")){
        pdf(paste0("figures/violin-",climate,"-",quant,".pdf"))
        .violin.plot(climate, quant, abs=FALSE)
        dev.off()
    }
}        

# Calculate the numbers
.percentage.stats <- function(index, var, clean, abs=FALSE, overshoot=FALSE, process=TRUE){
    index <- .simplify.quantiles(index, var, clean, abs)
    if(overshoot){
        index <- index > 1
    } else index <- index < 1 & index > 0
    if(process){
        percentages <- table(rowSums(index)) / nrow(index)
        percentages <- percentages[order(names(percentages))]
        something <- 1-percentages[1]
        at.least.two <- 1 - sum(percentages[1:2])
        at.least.three <- 1 - sum(percentages[1:3])
        everything <- percentages[6]
        return(setNames(
            c(something,at.least.two,at.least.three, everything),
            c("something","at.least.two","at.least.three","everything")))
    } else return(index)
}
.percentage.stats("bootstrap.index", "cld")
.percentage.stats("bootstrap.index", "frs")
.percentage.stats("bootstrap.index", "pet")
.percentage.stats("bootstrap.index", "pre")
.percentage.stats("bootstrap.index", "tmn")
.percentage.stats("bootstrap.index", "tmp")
.percentage.stats("bootstrap.index", "tmx")
.percentage.stats("bootstrap.index", "vap")
.percentage.stats("bootstrap.index", "wet")

combined <- .percentage.stats("bootstrap.index", "wet", process=FALSE) |
    .percentage.stats("bootstrap.index", "frs", process=FALSE) |
    .percentage.stats("bootstrap.index", "pet", process=FALSE) |
    .percentage.stats("bootstrap.index", "pre", process=FALSE) |
    .percentage.stats("bootstrap.index", "tmn", process=FALSE) |
    .percentage.stats("bootstrap.index", "tmp", process=FALSE) |
    .percentage.stats("bootstrap.index", "tmx", process=FALSE) |
    .percentage.stats("bootstrap.index", "vap", process=FALSE) |
    .percentage.stats("bootstrap.index", "wet", process=FALSE)
combined <- rowSums(combined)
combined <- table(combined) / length(combined)
1 - combined["0"]

.percentage.stats("bootstrap.index", "cld", overshoot=TRUE)
.percentage.stats("bootstrap.index", "frs", overshoot=TRUE)
.percentage.stats("bootstrap.index", "pet", overshoot=TRUE)
.percentage.stats("bootstrap.index", "pre", overshoot=TRUE)
.percentage.stats("bootstrap.index", "tmn", overshoot=TRUE)
.percentage.stats("bootstrap.index", "tmp", overshoot=TRUE)
.percentage.stats("bootstrap.index", "tmx", overshoot=TRUE)
.percentage.stats("bootstrap.index", "vap", overshoot=TRUE)
.percentage.stats("bootstrap.index", "wet", overshoot=TRUE)

o.combined <- .percentage.stats("bootstrap.index", "wet", process=FALSE) |
    .percentage.stats("bootstrap.index", "frs", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "pet", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "pre", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "tmn", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "tmp", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "tmx", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "vap", overshoot=TRUE, process=FALSE) |
    .percentage.stats("bootstrap.index", "wet", overshoot=TRUE, process=FALSE)
o.combined <- rowSums(o.combined)
o.combined <- table(o.combined) / length(o.combined)
1 - o.combined["0"]

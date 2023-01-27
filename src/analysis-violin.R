# Headers
source("src/headers.R")
stem <- "clnbin-clnspc-100-TRUE"

# Do the plots
.outline.plot <- function(stem, climate, quant, clean=NA, se.thresh=5, abs=FALSE, index="b.track.index", sims=FALSE){
    # Load and match data
    c(data,metadata) %<-% .load.indices(stem)
    sim.index <- as.data.frame(.simplify(data, index, quant, "spp.year.shuffle", clean, abs))
    sim.index$taxon <- metadata$taxon[match(rownames(sim.index), metadata$species)]
    index <- as.data.frame(.simplify(data, index, quant, "observed", clean, abs))
    index$taxon <- metadata$taxon[match(rownames(index), metadata$species)]
    data <- data.frame(index=index[,climate], taxon=index$taxon)
    sim.data <- data.frame(index=sim.index[,climate], taxon=sim.index$taxon)
    
    # Prepare colours and order for plotting
    cols <- setNames(c("forestgreen","forestgreen","red", "red", "red", "red","red"), c("plants", "fungi", "insects", "mammals", "birds", "amphibians", "reptiles"))
    data$taxon <- factor(data$taxon, levels=unique(as.character(data$taxon)))
    data$j.taxon <- as.numeric(factor(data$taxon))
    sim.data$taxon <- factor(sim.data$taxon, levels=unique(as.character(sim.data$taxon)))
    sim.data$j.taxon <- as.numeric(factor(sim.data$taxon))
    
    # Check for order of names as otherwise plots will be wrong
    if(!identical(names(cols), as.character(unique(data$taxon))))
        stop("Names of taxa in wrong order; plots will be wrong!")
    
    
    # Lookups for nicer names
    c.var <- setNames(c("cloud-cover","frost-days","potential-evapotranspiration","precipitation","min-temperature","mean-temperature","max-temperature","vapour-pressure"),
                      c("cld","frs","pet","pre","tmn","tmp","tmx","vap")
                      )[climate]
    q.var <- setNames(c("5th","25th","50th","75th","95th"), c("0.05","0.25","0.5","0.75","0.95"))[quant]

    # Calculate plot-points
    medians <- with(data, tapply(index, taxon, quantile, probs=.5))
    lower <- with(data, tapply(index, taxon, quantile, probs=.25))
    upper <- with(data, tapply(index, taxon, quantile, probs=.75))
    n.real <- with(data, tapply(index, taxon, length))
    s.medians <- with(sim.data, tapply(index, taxon, quantile, probs=.5))
    s.lower <- with(sim.data, tapply(index, taxon, quantile, probs=.25))
    s.upper <- with(sim.data, tapply(index, taxon, quantile, probs=.75))
    n.sim <- with(sim.data, tapply(index, taxon, length))

    
    par(mar=c(6.1,5.1,1.1,1.1))
    if(sims){
        xlims <- c(1,10.5)} else xlims <- c(1,7.5)
    eval(substitute(
        with(data, plot(index ~ jitter(j.taxon, factor=1.25), pch=20, col=alpha(cols,.25)[taxon], cex=.75, cex.lab=1.5, axes=FALSE, xlab="", xlim=xlims, ylab=expression(track[XXX](YYY)))), list(XXX=q.var, YYY=c.var)
    ))
    if(sims)
        with(sim.data, points(index ~ I(jitter(j.a.cat)+7), pch=20, col=alpha("blue",.25), cex=.25))
    
    arrows(1:7+.5, y0=lower, y1=upper, col=cols, length=0, lwd=5)
    arrows(1:7+.45, y0=medians, x1=1:7+.55, col=cols, length=0, lwd=10)
    arrows(1:7+.4, y0=s.lower, y1=s.upper, col="blue", length=0, lwd=5)
    arrows(1:7+.35, y0=s.medians, x1=1:7+.45, col="blue", length=0, lwd=10)
        
    abline(h=0, col="grey40", lwd=3, lty=2)
    abline(h=1, col="grey40", lwd=3, lty=2)
    
    # Do work
    axis(1, labels=names(cols), at=1:7, cex.axis=1.5)
    if(sims)
        axis(1, labels=unique(sim.data$alpha.cat), at=8:10, cex.axis=1.5)
    axis(2)
    if(sims){
        mtext(side=1, text=paste0("(n=",c(n.real,n.sim),")"), at=1:10, line=2)
    } else {
        mtext(side=1, text=paste0("(n=",n.real,")"), at=1:7, line=2)
    }

    # Add rasters
    if(sims){
        grid.raster(readPNG("phylopics/plant.png"),     .09, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/fungus.png"),    .18, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/insect.png"),    .27, .12, width=.02, just="bottom")
        grid.raster(readPNG("phylopics/mammal.png"),    .36, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/bird.png"),   .45, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/amphibian.png"), .54, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/reptile.png"),      .64, .12, width=.04, just="bottom")
    } else {
        grid.raster(readPNG("phylopics/plant.png"),     .095, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/fungus.png"),    .225, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/insect.png"),    .36, .12, width=.02, just="bottom")
        grid.raster(readPNG("phylopics/mammal.png"),    .49, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/bird.png"),   .62, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/amphibian.png"), .76, .12, width=.04, just="bottom")
        grid.raster(readPNG("phylopics/reptile.png"),      .89, .12, width=.04, just="bottom")
    }
    
}

for(climate in c("cld","frs","pet","pre","tmn","tmp","tmx","vap","wet")){
    for(quant in c("0.05","0.25","0.5","0.75","0.95")){
        pdf(paste0("figures/violin-",climate,"-",gsub(".","-",quant,fixed=TRUE),".pdf"), width=17, height=10)
        .outline.plot(stem, climate, quant, clean=5, index="b.track.index")
        dev.off()
    }
}        
pdf(paste0("figures/violin-tmp-0-5.pdf"), width=17, height=10);.outline.plot("tmp", "0.5", clean=5, sims=TRUE);dev.off()


# Calculate the numbers
.percentage.stats <- function(stem, index, var, clean, abs=FALSE, overshoot=FALSE, process=TRUE){
    browser()
    index <- .simplify.group(stem, index, null="observed", clean, abs)[,var]
    
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
.percentage.stats(stem, "b.track.index", "cld", clean=100)
.percentage.stats(stem, "b.track.index", "frs", clean=100)
.percentage.stats(stem, "b.track.index", "pet", clean=100)
.percentage.stats(stem, "b.track.index", "pre", clean=100)
.percentage.stats(stem, "b.track.index", "tmn", clean=100)
.percentage.stats(stem, "b.track.index", "tmp", clean=100)
.percentage.stats(stem, "b.track.index", "tmx", clean=100)
.percentage.stats(stem, "b.track.index", "vap", clean=100)
.percentage.stats(stem, "b.track.index", "wet", clean=100)

combined <- .percentage.stats(stem, "b.track.index", "wet", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "frs", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "pet", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "pre", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "tmn", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "tmp", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "tmx", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "vap", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "wet", process=FALSE, clean=100)
combined <- rowSums(combined)
combined <- table(combined) / length(combined)
1 - combined["0"]

.percentage.stats(stem, "b.track.index", "cld", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "frs", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "pet", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "pre", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "tmn", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "tmp", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "tmx", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "vap", overshoot=TRUE, clean=100)
.percentage.stats(stem, "b.track.index", "wet", overshoot=TRUE, clean=100)

o.combined <- .percentage.stats(stem, "b.track.index", "wet", process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "frs", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "pet", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "pre", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "tmn", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "tmp", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "tmx", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "vap", overshoot=TRUE, process=FALSE, clean=100) |
    .percentage.stats(stem, "b.track.index", "wet", overshoot=TRUE, process=FALSE, clean=100)
o.combined <- rowSums(o.combined)
o.combined <- table(o.combined) / length(o.combined)
1 - o.combined["0"]

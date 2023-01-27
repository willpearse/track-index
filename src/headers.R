####################################
# Packages #########################
####################################
if(!require("data.table")){
    install.packages("data.table"); library(data.table)}
if(!require("raster")){
    install.packages("raster"); library(raster)}
if(!require("RNetCDF")){
    install.packages("RNetCDF"); library(RNetCDF)}
if(!require("parallel")){
    install.packages("parallel"); library(parallel)}
if(!require("phytools")){
    install.packages("phytools"); library(phytools)}
if(!require("caper")){
    install.packages("caper"); library(caper)}
if(!require("abind")){
    install.packages("abind"); library(abind)}
if(!require("ks")){
    install.packages("ks"); library(ks)}
if(!require("geiger")){
    install.packages("geiger"); library(geiger)}
if(!require("suppdata")){
    install.packages("suppdata"); library(suppdata)}
if(!require("pez")){
    install.packages("pez"); library(pez)}
if(!require("scales")){
    install.packages("scales"); library(scales)}
if(!require("reldist")){
    install.packages("reldist"); library(reldist)}
if(!require("gdata")){
    install.packages("gdata"); library(gdata)}
if(!require("corrplot")){
    install.packages("corrplot"); library(corrplot)}
if(!require("viridis")){
    install.packages("viridis"); library(viridis)}
if(!require("plotrix")){
    install.packages("plotrix"); library(plotrix)}
if(!require("png")){
    install.packages("png"); library(png)}
if(!require("grid")){
    install.packages("grid"); library(grid)}
if(!require("zeallot")){
    install.packages("zeallot"); library(zeallot)}

####################################
# Global settings ##################
####################################
gbif.dir <- "/home/will/bear-share/gbif/"
if(gbif.dir!="/home/will/bear-share/gbif/")
    warning("GBIF data is not in default directory: ensure `download-gbif.rb` is similarly updated")

cru.dir <- "/home/will/bear-share/cru/"
if(cru.dir!="/home/will/bear-share/cru/")
    warning("CRU data is not in default directory: ensure `download-cru.rb` is similarly updated")

output.dir <- "clean-data/"

quantiles <- c(.05, .25, .5, .75, .95)

options(mc.cores=24)
if(options("mc.cores")==24)
    warning("Assuming you have 24 cores on your computer: if not, change `headers.R`")

data.options <- as.data.frame(expand.grid(
    min.records=c(500, 1000),
    clean.binomial=c(TRUE,FALSE),
    clean.spatial=c(TRUE,FALSE),
    binary=c(TRUE, FALSE)
))

data.options$output.dir <- with(data.options,
                                paste0("clean-data/",
                                       paste(
                                           ifelse(clean.binomial, "clnbin", "rawbin"),
                                           ifelse(clean.spatial, "clnspc", "rawspc"),
                                           min.records,
                                           binary,
                                           sep="-"
                                       )), "/")

####################################
# Functions ########################
####################################
prog.bar <- function(x, y){
    if(y < 100){
        cat(".")} else {
            z <- Filter(function(z) z>=0, seq(1,y,length.out=100)-x)
            if(length(z) > 0)
                tryCatch(if(z[1] < 1) if((length(z) %% 10)==0) cat("|") else cat("."), error=function(z) cat("."))
        }
}    
.load.gbif <- function(obs, specimens, min.records=2000, clean.binomial=FALSE, clean.spatial=FALSE){
    data <- data.table(scientificName="", decimalLatitude=0.0, decimalLongitude=0.0, year=0, type="")
    data$issue <- ""
    if(!missing(obs)){
        tmp <- fread(obs)
        tmp$type <- "observations"
        data <- rbind(data, tmp)
    }
    if(!missing(specimens)){
        tmp <- fread(specimens)
        tmp$type <- "specimens"
        data <- rbind(data, tmp)
    }
    if(missing(specimens) & missing(obs))
        stop("Must supply *some* data to load!...")
    data <- data[-1,]

    if(clean.spatial)
        data <- data[issue == "",]
    
    if(clean.binomial){
        spp <- unique(data$scientificName)
        spp <- setNames(tolower(sapply(strsplit(spp, " "), function(x) paste(x[1:2], collapse="_"))), spp)
        data$scientificName <- spp[data$scientificName]
    }

    # Find most common species and subset
    spp.table <- sort(table(data$scientificName), TRUE)
    spp.table <- names(spp.table)[spp.table >= min.records]
    data <- data[scientificName %in% spp.table,]

    # Subset by year
    data <- data[year > 1955 & year <= 2015,]
    return(data)
}
.assign.outcome <- function(present, past, projected){
    if(present==past)
        return("track")
    if(present==projected)
        return("static-distribution")

    # Warming
    if(projected > past){
        if(present < projected){
            if(present > past)
                return("lag")
            return("overshoot")
        }
        if(present > projected){
            return("undershoot")
        }
    }
    
    # Cooling
    if(projected < past){
        if(present > projected){
            if(present < past)
                return("lag")
            return("overshoot")
        }
        if(present < projected){
            return("undershoot")
        }
    }
    return("ERROR")
}
.simplify <- function(data, index, quantile, null, clean=NA, abs=FALSE){
    data <- t(sapply(data, function(x) x[null,,quantile,index]))
    if(!is.na(clean)){
        keep <- rep(TRUE, nrow(data))
        for(i in seq_len(ncol(data)))
            keep <- keep & abs(data[,i]) <= clean & (!is.na(data[,i]))
        data <- data[keep,]
    }
    if(abs)
        data <- abs(data)
    return(data)
}
.simplify.group <- function(stem, index, quantile, null="observed", clean=NA, abs=FALSE)
    return(rbind(
        data.frame(.simplify(.load.taxon(stem,"plants"), index, quantile, null, clean, abs), taxon="plants"),
        data.frame(.simplify(.load.taxon(stem,"fungi"), index, quantile, null, clean, abs), taxon="fungi"),
        data.frame(.simplify(.load.taxon(stem,"insects"), index, quantile, null, clean, abs), taxon="insects"),
        data.frame(.simplify(.load.taxon(stem,"mammals"), index, quantile, null, clean, abs), taxon="mammals"),
        data.frame(.simplify(.load.taxon(stem,"birds"), index, quantile, null, clean, abs), taxon="birds"),
        data.frame(.simplify(.load.taxon(stem,"reptiles"), index, quantile, null, clean, abs), taxon="reptiles"),
        data.frame(.simplify(.load.taxon(stem,"amphibians"), index, quantile, null, clean, abs), taxon="amphibians")
    ))
.simplify.quantile <- function(stem, index, var, clean=NA, abs=FALSE)
    return(rbind(
        data.frame(.simplify(.load.taxon(stem,"plants"), index, quantile, clean, abs), taxon="plants"),
        data.frame(.simplify(.load.taxon(stem,"fungi"), index, quantile, clean, abs), taxon="fungi"),
        data.frame(.simplify(.load.taxon(stem,"insects"), index, quantile, clean, abs), taxon="insects"),
        data.frame(.simplify(.load.taxon(stem,"mammals"), index, quantile, clean, abs), taxon="mammals"),
        data.frame(.simplify(.load.taxon(stem,"birds"), index, quantile, clean, abs), taxon="birds"),
        data.frame(.simplify(.load.taxon(stem,"reptiles"), index, quantile, clean, abs), taxon="reptiles"),
        data.frame(.simplify(.load.taxon(stem,"amphibians"), index, quantile, clean, abs), taxon="amphibians")
    ))

.load.taxon <- function(stem, x){
    output <- readRDS(paste0("clean-data/",stem,x,"-index.RDS"))
    output <- Filter(Negate(is.character), output)
    #output <- Filter(function(y) all(is.na(y)), output)
    return(output)
}
.load.indices <- function(stem="rawbin-clnspc-1000-TRUE"){
    plants <- .load.taxon(stem, "plants")
    fungi <- .load.taxon(stem, "fungi")
    insects <- .load.taxon(stem, "insects")
    mammals <- .load.taxon(stem, "mammals")
    reptiles <- .load.taxon(stem, "reptiles")
    birds <- .load.taxon(stem, "birds")
    amphibians <- .load.taxon(stem, "amphibians")
    combined <- unlist(list(plants, fungi, insects, mammals, birds, amphibians, reptiles), recursive=FALSE)
    array <- combined[sapply(combined, class) == "array"]
    metadata <- data.frame(
        species=names(array),
        taxon=rep(
            c("plants", "fungi", "insects", "mammals", "birds", "amphibians", "reptiles"),
            sapply(list(plants,fungi,insects,mammals,birds,amphibians,reptiles), length)
        )
    )
    return(list(index=array, metadata=metadata))
}

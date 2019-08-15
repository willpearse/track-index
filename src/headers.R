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

####################################
# Global settings ##################
####################################
gbif.dir <- "~/bear-share/gbif/"
if(gbif.dir!="~/bear-share/gbif/")
    warning("GBIF data is not in default directory: ensure `download-gbif.rb` is similarly updated")

cru.dir <- "~/bear-share/cru/"
if(gbif.dir!="~/bear-share/cru/")
    warning("CRU data is not in default directory: ensure `download-cru.rb` is similarly updated")

output.dir <- "clean-data/"

options(mc.cores=24)
if(options("mc.cores"==24))
    warning("Assuming you have 24 cores on your computer: if not, change `headers.R`")

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
.load.gbif <- function(obs, specimens, min.records=2000){
    data <- data.table(scientificname="", decimallatitude=0.0, decimallongitude=0.0, year=0, type="")
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

    # Find most common species and subset
    spp.table <- sort(table(data$scientificname), TRUE)
    spp.table <- names(spp.table)[spp.table >= min.records]
    data <- data[scientificname %in% spp.table,]

    # Subset by year
    data <- data[year > 1955 & year <= 2015,]
    return(data)
}
.create.raster.stack <- function(x, long=seq(-179.75,179.75,by=0.5), lat=seq(-89.75,89.75,by=0.5)){
    output <- raster(t(x[,,1][,ncol(x):1]), xmn=min(long), xmx=max(long), ymn=min(lat), ymx=max(lat))
    base <- raster(xmn=min(long), xmx=max(long), ymn=min(lat), ymx=max(lat))
    for(i in seq(2, dim(x)[3]))
        output <- stack(output, raster(t(x[,,i][,ncol(x):1]), xmn=min(long), xmx=max(long), ymn=min(lat), ymx=max(lat)))
    return(output)
}
.calc.track <- function(i){
    name <- paste0("taxon_", i, "_results.RDS")
    data <- wrap.list[[i]]
    quantiles <- seq(.1,.9,.05)
    results <- array(dim=c(length(unique(data$scientificname)),ncol=14,length(quantiles)), dimnames=list(unique(data$scientificname), c("present","present.mad","past","past.mad","projected","projected.mad","n.pre","n.post","r.pre","r.post","r.proj","r.total","pre.lat","post.lat"), quantiles))
    for(i in seq_len(nrow(results))){
        # Pre-allocation and calculation
        curr.data <- data[scientificname == rownames(results)[i],]
        post.data <- curr.data[year %in% 1990:2015,]
        pre.data <- curr.data[year %in% 1955:1980,]
        post.years <- unique(post.data$year)
        post.data$env <- -999
        for(k in seq_along(post.years))
            post.data[year==post.years[k],]$env <- extract(temp[[post.data$year[k]-1900]], post.data[year==post.years[k],.(decimallongitude,decimallatitude)])

        # Chuck out if potential sampling problems
        if(nrow(pre.data) == 0 | nrow(post.data) == 0)
            next
        
        # Past and projected data extraction
        pre.data$pre.env <- -999
        pre.data$post.env <- -999
        pre.years <- unique(pre.data$year)
        for(k in seq_along(pre.years)){
            pre.data[year==pre.years[k],]$pre.env <- extract(temp[[pre.data$year[k]-1900]], pre.data[year==pre.years[k],.(decimallongitude,decimallatitude)])
            pre.data[year==pre.years[k],]$post.env <- extract(temp[[pre.data$year[k]-1900+25]], pre.data[year==pre.years[k],.(decimallongitude,decimallatitude)])
        }

        # Calculate differently depending on sampling
        if(nrow(pre.data) < nrow(post.data)){
            # More data post-1980 than before; sub-sample accordingly
            results[i,"past",] <- quantile(pre.data$pre.env, probs=quantiles, na.rm=TRUE)
            results[i,"projected",] <- quantile(pre.data$post.env, probs=quantiles, na.rm=TRUE)
            results[i,"past.mad",] <- results[i,"projected.mad",] <- 0
            t <- replicate(100, sample(post.data$env, nrow(pre.data)), simplify=FALSE)
            results[i,"present",] <- apply(sapply(t, quantile, probs=quantiles, na.rm=TRUE), 1, median)
            results[i,"present.mad",] <- apply(sapply(t, quantile, probs=quantiles, na.rm=TRUE), 1, mad)
        } else {
            if(nrow(pre.data) > nrow(post.data)){
                # More data pre-1980 than post; sub-sample accordingly
                t <- replicate(100, sample(pre.data$pre.env, nrow(post.data)), simplify=FALSE)
                results[i,"past",] <- apply(sapply(t, quantile, probs=quantiles, na.rm=TRUE), 1, median)
                results[i,"past.mad",] <- apply(sapply(t, quantile, probs=quantiles, na.rm=TRUE), 1, mad)
                t <- replicate(100, sample(pre.data$post.env, nrow(post.data)), simplify=FALSE)
                results[i,"projected",] <- apply(sapply(t, quantile, probs=quantiles, na.rm=TRUE), 1, median)
                results[i,"projected.mad",] <- apply(sapply(t, quantile, probs=quantiles, na.rm=TRUE), 1, sd)
                results[i,"present",] <- quantile(post.data$env, probs=quantiles, na.rm=TRUE)
                results[i,"present.mad",] <- 0
            } else {
                # What are the odds? The same amount of data in each period
                results[i,"past",] <- quantile(pre.data$pre.env, probs=quantiles, na.rm=TRUE)
                results[i,"projected",] <- quantile(pre.data$post.env, probs=quantiles, na.rm=TRUE)
                results[i,"present",] <- quantile(post.data$env, probs=quantiles, na.rm=TRUE)
                results[i,"past.mad",] <- results[i,"projected.mad",] <- results[i,"present.mad",] <- 0
            }
        }

        results[i,"n.pre",] <- nrow(pre.data)
        results[i,"n.post",] <- nrow(post.data)

        results[i,"r.pre",] <- diff(range(pre.data$pre.env, na.rm=TRUE))
        results[i,"r.proj",] <- diff(range(pre.data$post.env, na.rm=TRUE))
        results[i,"r.post",] <- diff(range(post.data$env, na.rm=TRUE))
        results[i,"r.total",] <- diff(range(c(post.data$env,pre.data$post.env,pre.data$pre.env), na.rm=TRUE))

        results[i,"pre.lat",] <- median(pre.data$decimallatitude, na.rm=TRUE)
        results[i,"post.lat",] <- median(post.data$decimallatitude, na.rm=TRUE)
    }
    saveRDS(results, name)
    return(results)
}

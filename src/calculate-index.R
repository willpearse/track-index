# Headers
source("src/headers.R")

# Link past and present distributions with past and present climate
map.climate <- function(gbif.data, env, pres.years, past.years, yr.diff, binary){
    .get <- function(...){
        if(binary){
            cells <- mask(env, rasterize(points, env, fun="first"))@data@values
            return(cells[!is.na(cells)])
        }
        points <- extract(env, gbif.data[,.(decimallongitude,decimallatitude)])@data@values
        return(points[!is.na(points)])
    }
    # Pre-allocate
    pres.dis.pres.env <- pres.dis.past.env <- vector("list", length(pres.years)) 
    past.dis.pres.env <- past.dis.past.env <- vector("list", length(past.years))

    # Extract
    for(i in seq_along(pres.years)){
        pres.dis.pres.env[[i]] <- .get(gbif.data[year==pres.years[i]], env[[pres.years[i]-1900]])
        pres.dis.past.env[[i]] <- .get(gbif.data[year==pres.years[i]], env[[pres.years[i]-1900-yr.diff]])
    }
        
    for(i in seq_along(past.years)){
        past.dis.pres.env[[i]] <- .get(gbif.data[year==past.years[i]], env[[past.years[i]-1900+yr.diff]])
        past.dis.past.env[[i]] <- .get(gbif.data[year==past.years[i]], env[[past.years[i]-1900]])
    }
    
    # Reformat and return
    return(list(
        pres.dis.pres.env=unlist(pres.dis.pres.env), pres.dis.past.env=unlist(pres.dis.past.env),
        past.dis.pres.env=unlist(past.dis.pres.env), past.dis.past.env=unlist(past.dis.past.env)
    ))
}

# Calculate index
calc.metric <- function(pres.dis.pres.env, pres.dis.past.env, past.dis.pres.env, past.dis.past.env, quantile=.5, n.boot=999){
    # Calculate observed quantiles and indices
    quants <- sapply(list(pres.dis.pres.env, pres.dis.past.env, past.dis.pres.env, past.dis.past.env), quantile, probs=quantile, na.rm=TRUE)
    track.index <- (quants["pres.dis.pres.env"]-quants["past.dis.past.env"]) / (quants["past.dis.pres.env"]-quants["past.dis.past.env"])
    null.index <- (quants["pres.dis.pres.env"]-quants["pres.dis.past.env"]) / (quants["past.dis.pres.env"]-quants["past.dis.past.env"])

    # Bootstrap resampling
    .bootstrap <- function(x, n, n.boot, quantile){
        x <- replicate(n.boot, sample(x, n, replace=TRUE), simplify=FALSE)
        return(sapply(x, quantile, probs=quantile, na.rm=TRUE))
    }
    b.pres.dis.pres.env <- .bootstrap(pres.dis.pres.env, length(past.dis.past.env), n.boot, quantile)
    b.pres.dis.past.env <- .bootstrap(pres.dis.past.env, length(past.dis.past.env), n.boot, quantile)
    b.past.dis.pres.env <- .bootstrap(past.dis.pres.env, length(pres.dis.pres.env), n.boot, quantile)
    b.past.dis.past.env <- .bootstrap(past.dis.past.env, length(pres.dis.pres.env), n.boot, quantile)

    # Calculate bootstrapped indices and summaries
    b.track.index <- (b.pres.dis.pres.env-b.past.dis.past.env) / (b.past.dis.pres.env-b.past.dis.past.env)
    mad.track.index <- mad(b.track.index)
    b.track.index <- median(b.track.index)
    b.null.index <- (b.pres.dis.pres.env-b.past.dis.past.env) / (b.pres.dis.pres.env-b.past.dis.past.env)
    mad.null.index <- mad(b.null.index)
    b.null.index <- median(b.null.index)
    
    # Format and return
    return(c(
        track.index=track.index, b.track.index=b.track.index, mad.track.index=mad.track.index,
        null.index=null.index, b.null.index=b.null.index, mad.null.index=mad.null.index,
        pres.dis.pres.env=quants["pres.dis.pres.env"], pres.dis.past.env=quants["pres.dis.past.env"], past.dis.pres.env=quants["past.dis.pres.env"], past.dis.past.env=quants["past.dis.past.env"]
    ))
}

# Main worker for observed and null data
calc.metric.parallel <- function(data, quantiles, rasters, taxon, output.dir, binary, null=c("obs","shift","niche","sampling")){
    .int.calc <- function(sp, binary){
        # Setup and subset data
        curr.data <- data[scientificname == sp,]
        output <- array(
            dim=c(2, length(rasters), length(quantiles), 7),
            dimnames=list(c("observed","spp.year.shuffle")names(rasters), as.character(quantiles), c(track.index, b.track.index, mad.track.index, null.index, b.null.index, mad.null.index, pres.dis.pres.env, pres.dis.past.env, past.dis.pres.env, past.dis.past.env))
        )
        
        # Do work per raster
        for(i in seq_along(rasters)){
            t <- map.climate(curr.data, rasters[[i]], 1990:2015, 1955:1980, 35, binary)
            output[1,i,,] <- t(sapply(quantiles, function(x) calc.metric(t$pres.dis.pres.env, t$pres.dis.past.env, t$past.dis.pres.env, t$past.dis.past.env, x)))
            curr.data$year <- sample(curr.data$year)
            t <-.map.climate(curr.data, rasters[[i]], 1990:2015, 1955:1980, 35, binary)
            output[2,i,,] <- t(sapply(quantiles, function(x) calc.metric(t$pres.dis.pres.env, t$pres.dis.past.env, t$past.dis.pres.env, t$past.dis.past.env, x)))
        }
        return(output)
    }
    
    return(mcMap(.int.calc, unique(data$scientificname, binary=binary)))
}

do.work <- function(i){
    .taxon <- function(file, taxon){
        gbif <- .load.gbif(obs=paste0(gbif.dir,file), min.records=opts["min.records"], clean.binomial=opts["clean.binomial"], clean.gbif=opts["clean.gbif"])
        saveRDS(
            .calc.track.parallel(gbif, quantiles, cru, taxon, opts["output.dir"], binary=opts["binary"]),
            paste0(opts["output.dir"],paste0(taxon,"-index.RDS"))
        )
    }

    opts <- data.options[i,]
    .taxon("amphibia/0011411-200221144449610_abbreviated.txt", "amphibians")
    .taxon("asco-basidio/0011412-200221144449610_abbreviated.txt", "fungi")
    .taxon("aves/0011413-200221144449610_abbreviated.txt", "birds")
    .taxon("insecta/0012312-200221144449610_abbreviated.txt", "insects")
    .taxon("mammalia/0012314-200221144449610_abbreviated.txt", "mammals")
    .taxon("plantae/0012316-200221144449610_abbreviated.txt", "plants")
    .taxon("reptillia/0012360-200221144449610_abbreviated.txt", "reptiles")
}


# Load climate rasters
cru <- list(
    cld=readRDS(paste0(output.dir,"cru-cld.RDS")),
    frs=readRDS(paste0(output.dir,"cru-frs.RDS")),
    pet=readRDS(paste0(output.dir,"cru-pet.RDS")),
    pre=readRDS(paste0(output.dir,"cru-pre.RDS")),
    tmn=readRDS(paste0(output.dir,"cru-tmn.RDS")),
    tmp=readRDS(paste0(output.dir,"cru-tmp.RDS")),
    tmx=readRDS(paste0(output.dir,"cru-tmx.RDS")),
    vap=readRDS(paste0(output.dir,"cru-vap.RDS")),
    wet=readRDS(paste0(output.dir,"cru-wet.RDS"))
)    

# Do work across all data formatting types
for(i in seq_len(nrow(data.options)))
    do.work(i)

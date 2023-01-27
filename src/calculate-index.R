# Headers
source("src/headers.R")

# Link past and present distributions with past and present climate
map.climate <- function(gbif.data, env, pres.years, past.years, yr.diff, binary, n.boots){
  .get <- function(d, e){
    if(nrow(d)==0)
      return(numeric(0))
    if(binary){
      cells <- mask(e, rasterize(d[,.(decimalLongitude,decimalLatitude)], e, fun="first"))@data@values
      return(cells[!is.na(cells)])
    }
    points <- extract(e, d[,.(decimalLongitude,decimalLatitude)])
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
  if(length(pres.dis.pres.env)==0 | length(pres.dis.past.env)==0 | length(past.dis.pres.env)==0 | length(past.dis.past.env)==0){
    return(c(
      track.index=NA, stat.index=NA, 
      pres.dis.pres.env=NA, pres.dis.past.env=NA, past.dis.pres.env=NA, past.dis.past.env=NA
    ))
  }
  
  # Calculate observed quantiles and indices
  quants <- setNames(sapply(list(pres.dis.pres.env, pres.dis.past.env, past.dis.pres.env, past.dis.past.env), quantile, probs=quantile, na.rm=TRUE), c("pres.dis.pres.env", "pres.dis.past.env", "past.dis.pres.env", "past.dis.past.env"))
  track.index <- (quants["pres.dis.pres.env"]-quants["past.dis.past.env"]) / (quants["past.dis.pres.env"]-quants["past.dis.past.env"])
  stat.index <- (quants["past.dis.pres.env"]-quants["pres.dis.pres.env"]) / (quants["past.dis.pres.env"]-quants["past.dis.past.env"])

  # Format and return
  return(c(
    track.index=track.index, stat.index=stat.index, pres.dis.pres.env=quants["pres.dis.pres.env"], pres.dis.past.env=quants["pres.dis.past.env"], past.dis.pres.env=quants["past.dis.pres.env"], past.dis.past.env=quants["past.dis.past.env"]
  ))
}

# Main worker for observed and null data
calc.metric.parallel <- function(data, quantiles, rasters, taxon, output.dir, binary){
  .int.calc <- function(sp, binary){
    # Setup and subset data
    curr.data <- data[scientificName == sp,]
    output <- array(
      dim=c(2, length(rasters), length(quantiles), 10),
      dimnames=list(c("observed","stat.null","track.nulk"), names(rasters), as.character(quantiles), c("track.index", "stat.index", "pres.dis.pres.env", "pres.dis.past.env", "past.dis.pres.env", "past.dis.past.env"))
    )
    
    # Do work per raster
    for(i in seq_along(rasters)){
      # Observed metric
      t <- map.climate(curr.data, rasters[[i]], 1990:2015, 1955:1980, 35, binary)
      output[1,i,,] <- t(sapply(quantiles, function(x) calc.metric(t$pres.dis.pres.env, t$pres.dis.past.env, t$past.dis.pres.env, t$past.dis.past.env, x)))

      # Stationary null
      s.data <- curr.data
      holder <- matrix(NA, nrow=n.boot, ncol=length(quantiles))
      merged.pres.past <- c(t$pres.dis.pres.env, t$past.dis.past.env)
      for(j in seq_len(n.boot)){
        s.data$year <- sample(s.data$year, size=nrow(s.data))
        t <- map.climate(s.data, rasters[[i]], 1990:2015, 1955:1980, 35, binary)
        holder[j,] <- t(sapply(quantiles, function(x) calc.metric(t$pres.dis.pres.env, t$pres.dis.past.env, t$past.dis.pres.env, t$past.dis.past.env, x)))[,"track.index"]
      }
      for(j in seq_len(ncol(holder)))
        output[2,i,1,j] <- rank(c(holder[,j], output[1,i,,"track.index"]))

      # Tracking null
      s.data <- curr.data
      t <- map.climate(s.data, rasters[[i]], 1990:2015, 1955:1980, 35, binary)
      holder <- matrix(NA, nrow=n.boot, ncol=length(quantiles))
      merged.pres.past <- c(t$pres.dis.pres.env, t$past.dis.past.env)
      for(j in seq_len(n.boot)){
        merged.pres.past <- sample(merged.pres.past, replace=FALSE)
        holder[i,] <- t(sapply(quantiles, function(x) calc.metric(merged.pres.past[seq(1,length(t$pres.dis.pres.env))], t$pres.dis.past.env, t$past.dis.pres.env, merged.pres.past[seq(length(t$pres.dis.pres.env)+1,length(merged.pres.past))], x)))[,"track.index"]
      }
      for(j in seq_len(ncol(holder)))
        output[3,i,1,j] <- rank(c(holder[,j], output[1,i,,"track.index"]))
    }
    return(output)
  }
  return(mcMap(.int.calc, unique(data$scientificName), binary=binary))
}

do.work <- function(i){
  .taxon <- function(file, taxon){
    gbif <- .load.gbif(obs=paste0(gbif.dir,file), min.records=unlist(opts["min.records"]), clean.binomial=unlist(opts["clean.binomial"]), clean.spatial=unlist(opts["clean.spatial"]))
    saveRDS(
      calc.metric.parallel(gbif, quantiles, cru, taxon, unlist(opts["output.dir"]), binary=unlist(opts["binary"])),
      paste0(unlist(opts["output.dir"]),paste0(taxon,"-index.RDS"))
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

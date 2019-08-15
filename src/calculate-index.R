# Headers
source("src/headers.R")

# Load and clean GBIF data
plants <- .load.gbif(
    obs=paste0(gbif.dir,"plantae/0001849-180131172636756_abbreviated.txt"),
    specimens=paste0(gbif.dir,"plantae/0001851-180131172636756_abbreviated.txt")
)
mammals <- .load.gbif(obs=paste0(gbif.dir,"mammalia/0000171-180412121330197_abbreviated.txt"))
amphibians <- .load.gbif(obs=paste0(gbif.dir,"amphibia/0034734-180508205500799_abbreviated.txt"))
fungi <- .load.gbif(obs=paste0(gbif.dir,"asco-basidio/0034724-180508205500799_abbreviated.txt"))
birds <- .load.gbif(obs=paste0(gbif.dir,"aves/0000172-180412121330197_abbreviated.txt"))
insects <- .load.gbif(obs=paste0(gbif.dir,"insecta/0034726-180508205500799_abbreviated.txt"))
reptiles <- .load.gbif(obs=paste0(gbif.dir,"reptillia/0034728-180508205500799_abbreviated.txt"))

# Load and clean climate data
temp <- readRDS(paste0(output.dir,"cru-temp.RDS"))

# Do work and save results
wrap.list <- list(plants, mammals, amphibians, fungi, birds, insects, reptiles)
results <- mcMap(.calc.track, seq_along(wrap.list))
saveRDS(results, paste0(output.dir,"cru-temp.RDS"))

# Headers
source("src/headers.R")

# Wrapper functions
.univariate <- function(data, explanatory){
    tests <- lapply(data, function(x) cor.test(x, explanatory))
    cors <- sapply(tests, function(x) x$estimate)
    p <- sapply(tests, function(x) x$p.value)
    p <- ifelse(p < .2, round(p,3), NA)
    output <- matrix(c(cors,p), nrow=2, byrow=TRUE)
    colnames(output) <- names(data); rownames(output) <- c("r","p")
    return(output)
}

# Load data
plants <- readRDS("clean-data/plants-index.RDS")
fungi <- readRDS("clean-data/plants-index.RDS")
insects <- readRDS("clean-data/plants-index.RDS")
mammals <- readRDS("clean-data/mammals-index.RDS")
reptiles <- readRDS("clean-data/reptiles-index.RDS")
birds <- readRDS("clean-data/birds-index.RDS")

neotoma <- read.csv("http://www.esapubs.org/archive/ecol/E096/269/Data_Files/Amniote_Database_Aug_2015.csv", as.is=TRUE)

# Neotoma
neotoma$binomial <- with(neotoma, tolower(paste(genus, species, sep="_")))
combined.index$log.mass <- neotoma$adult_body_mass_g[match(rownames(combined.index), neotoma$binomial)]
combined.index$log.mass[combined.index$mass==-999] <- NA
combined.index$log.mass <- log(combined.index$log.mass)
.univariate(combined.index[,1:9], combined.index$log.mass)

#Glopnet
library(gdata)
glopnet <- read.xls("raw-data/glopnet.xls")
combined.index$species <- rownames(combined.index)
glopnet <- merge(glopnet, combined.index, by.x="Species", by.y="species")
.univariate(glopnet[,names(combined.index)[1:9]], glopnet$log.LL)
.univariate(glopnet[,names(combined.index)[1:9]], glopnet$log.LMA)
.univariate(glopnet[,names(combined.index)[1:9]], glopnet$log.Nmass)
.univariate(glopnet[,names(combined.index)[1:9]], glopnet$log.Amass)


elton.birds <- read.delim("raw-data/BirdFuncDat.txt", as.is=TRUE)
elton.birds$Scientific <- tolower(gsub(" ", "_", elton.birds$Scientific))
elton.mam <- read.delim("raw-data/MamFuncDat.txt", as.is=TRUE)
elton.mam$Scientific <- tolower(gsub(" ", "_", elton.mam$Scientific))
elton <- setNames(c(elton.mam$BodyMass.Value, elton.birds$BodyMass.Value), c(elton.mam$Scientific,elton.birds$Scientific))
combined.index$elton.log.mass <- log(elton[combined.index$species])
.univariate(combined.index[,1:9], combined.index$elton.log.mass)
#... this gives the only thing I think is close: pre ~ elton.log.mass

#... playing around a little with only the present-day observations,
#and I get something here for the body mass things. So I think this is
#reassuring because it means there's some kind of pattern...


# Simulating PCA stuff
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



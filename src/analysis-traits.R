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
.combine.data <- function(index, quant, clean=10, abs=FALSE){
    data <- .simplify.group(index, quant, clean, abs)
    data$species <- rownames(data)

    data$bm.neotoma <- neotoma$log.mass[match(data$species, neotoma$binomial)]
    data$bm.elton <- elton.bm[data$species]

    data$log.ll <- glopnet$log.LL[match(data$species, glopnet$Species)]
    data$log.lma <- glopnet$log.LMA[match(data$species, glopnet$Species)]
    data$log.Nmass <- glopnet$log.Nmass[match(data$species, glopnet$Species)]
    data$log.Amass <- glopnet$log.Amass[match(data$species, glopnet$Species)]
    return(data)    
}
.corr.mat <- function(index, quant, clean, abs=FALSE){
    data <- .combine.data(index, quant, clean, abs)

    c.var <- setNames(
        c("cloud-cover","frost-days","potential-evapotranspiration","precipitation","min-temperature","mean-temperature","max-temperature","vapour-pressure"),
                      c("cld","frs","pet","pre","tmn","tmp","tmx","vap")
    )
    q.var <- setNames(
        c("5th","25th","50th","75th","95th"),
        c("0.05","0.25","0.5","0.75","0.95")
    )[quant]
    t.var <- setNames(
        c("body-mass","body-mass","leaf-lifespan","leaf-mass/area","leaf-N","photosynthetic-capacity"),
        c("bm.neotoma","bm.elton","log.ll","log.lma","log.Nmass","log.Amass")
    )
    
    mat <- matrix(nrow=6, ncol=8, dimnames=list(t.var, c.var))
    for(i in seq_along(c.var)){
        for(j in seq_along(t.var)){
            mat[j,i] <- cor.test(data[,names(c.var)[i]], data[,names(t.var)[j]])$estimate
        }
    }
    return(mat)
}
.plot.corr.mat <- function(index, comparison, quant, clean, abs=FALSE){
    index.mat <- t(.corr.mat(index, quant, clean, abs))
    comparison.mat <- t(.corr.mat(comparison, quant))
    cols <- colorRampPalette(c("red", "white", "blue"))
                   
    comparison.cuts <- as.numeric(cut(as.numeric(comparison.mat), breaks=200))
    comparison.cuts <- cols(200)[comparison.cuts]
    index.cuts <- as.numeric(cut(as.numeric(index.mat), breaks=200))
    index.cuts <- cols(200)[index.cuts]
    dummy.mat <- matrix(0, nrow=nrow(index.mat), ncol=ncol(index.mat), dimnames=dimnames(index.mat))
    
    corrplot(comparison.mat, method="square", is.corr=FALSE, cl.lim=c(-1,1), tl.srt=45)
    corrplot(comparison.mat, method="square", is.corr=FALSE, cl.lim=c(-1,1), tl.srt=45, bg=alpha("white",0), add=TRUE, addgrid.col=NA, col=colorRampPalette(c("blue","white", "blue"))(200), cl.pos="n", tl.pos="n")
    corrplot(index.mat, method="square", is.corr=FALSE, cl.lim=c(-1,1), tl.srt=45, bg=alpha("white",0), col=colorRampPalette(c("red","white", "red"))(200), add=TRUE, addgrid.col=NA, cl.pos="n", , tl.pos="n")
}

# Load tracking data
plants <- readRDS("clean-data/plants-index.RDS")
fungi <- readRDS("clean-data/plants-index.RDS")
insects <- readRDS("clean-data/plants-index.RDS")
mammals <- readRDS("clean-data/mammals-index.RDS")
reptiles <- readRDS("clean-data/reptiles-index.RDS")
birds <- readRDS("clean-data/birds-index.RDS")
combined <- .simplify.group("bootstrap.index", "0.5", clean=10)

# Load Neotoma data
neotoma <- read.csv("raw-data/Amniote_Database_Aug_2015.csv", as.is=TRUE)
for(i in seq_along(names(neotoma)))
    neotoma[neotoma[,i]==-999,i] <- NA
neotoma$binomial <- with(neotoma, tolower(paste(genus, species, sep="_")))
neotoma$log.mass <- log10(neotoma$adult_body_mass_g)

# Load Glopnet
glopnet <- read.xls("raw-data/glopnet.xls")
glopnet$Species <- tolower(gsub(" ", "_", glopnet$Species))

# Load Elton
elton.birds <- read.delim("raw-data/BirdFuncDat.txt", as.is=TRUE)
elton.birds$Scientific <- tolower(gsub(" ", "_", elton.birds$Scientific))
elton.mam <- read.delim("raw-data/MamFuncDat.txt", as.is=TRUE)
elton.mam$Scientific <- tolower(gsub(" ", "_", elton.mam$Scientific))
elton.bm <- setNames(c(elton.mam$BodyMass.Value, elton.birds$BodyMass.Value), c(elton.mam$Scientific,elton.birds$Scientific))
elton.bm <- log10(elton.bm)

# Combine
pdf("figures/traits-05.pdf"); .plot.corr.mat("bootstrap.index", "present", "0.05", 5); dev.off()
pdf("figures/traits-25.pdf"); .plot.corr.mat("bootstrap.index", "present", "0.25", 5); dev.off()
pdf("figures/traits-50.pdf"); .plot.corr.mat("bootstrap.index", "present", "0.5", 5); dev.off()
pdf("figures/traits-75.pdf"); .plot.corr.mat("bootstrap.index", "present", "0.75", 5); dev.off()
pdf("figures/traits-95.pdf"); .plot.corr.mat("bootstrap.index", "present", "0.95", 5); dev.off()

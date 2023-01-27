# Headers
source("src/headers.R")
stem <- "clnbin-clnspc-100-TRUE"

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
.combine.data <- function(stem, index, quant, null, clean=10, abs=FALSE){
    c(data,metadata) %<-% .load.indices(stem)
    data <- as.data.frame(.simplify(data, index, quant, null, clean, abs))
    data$species <- rownames(data)
    data$taxon <- metadata$taxon[match(data$species, metadata$species)]
    
    data$bm.mamm <- neotoma$log.mass[match(data$species, neotoma$binomial)]
    data$bm.mamm[data$taxon != "mammals"] <- NA
    data$bm.bird <- elton.bm[data$species]
    data$bm.bird[data$taxon != "birds"] <- NA

    data$log.ll <- glopnet$log.LL[match(data$species, glopnet$Species)]
    data$log.lma <- glopnet$log.LMA[match(data$species, glopnet$Species)]
    data$log.Nmass <- glopnet$log.Nmass[match(data$species, glopnet$Species)]
    data$log.Amass <- glopnet$log.Amass[match(data$species, glopnet$Species)]
    return(data)    
}
.corr.mat <- function(stem, index, quant, null, clean, abs=FALSE, p.val=FALSE){
    data <- .combine.data(stem, index, quant, null, clean, abs)
    c.var <- setNames(
        c("clouds","frost","evapo-trans.","precipitation","min(temp)","mean(temp)","max(temp)","vapour","rainy day"),
        c("cld","frs","pet","pre","tmn","tmp","tmx","vap","wet")
    )
    q.var <- setNames(
        c("5th","25th","50th","75th","95th"),
        c("0.05","0.25","0.5","0.75","0.95")
    )[quant]
    t.var <- setNames(
        c("mammal mass","bird mass","leaf lifespan","leaf mass/area","leaf N","photosynth"),
        c("bm.mamm","bm.bird","log.ll","log.lma","log.Nmass","log.Amass")
    )
    
    mat <- matrix(nrow=length(t.var), ncol=length(c.var), dimnames=list(t.var, c.var))
    for(i in seq_along(c.var)){
        for(j in seq_along(t.var)){
            if(p.val){
                mat[j,i] <- cor.test(data[,names(c.var)[i]], data[,names(t.var)[j]])$p.value               
            } else {
                mat[j,i] <- cor.test(data[,names(c.var)[i]], data[,names(t.var)[j]])$estimate
            }
        }
    }
    return(mat)
}
.plot.corr.mat <- function(stem, index, comparison, quant, null, clean, abs=FALSE){
    index.mat <- t(.corr.mat(stem, index, quant, null, clean, abs))
    index.mat.p <- t(.corr.mat(stem, index, quant, null, clean, abs, p.val=TRUE))
    signif <- index.mat.p < .05
    index.mat.p <- matrix(0, nrow=nrow(index.mat), ncol=ncol(index.mat))
    index.mat.p[signif] <- 1
    comparison.mat <- t(.corr.mat(stem, comparison, quant, null, clean=NA))
    cols <- colorRampPalette(c("red", "white", "blue"))
                   
    comparison.cuts <- as.numeric(cut(as.numeric(comparison.mat), breaks=seq(-1,1,length.out=201)))
    comparison.cuts <- cols(200)[comparison.cuts]
    index.cuts <- as.numeric(cut(as.numeric(index.mat), breaks=seq(-1,1,length.out=201)))
    index.cuts <- cols(200)[index.cuts]
    dummy.mat <- matrix(0, nrow=nrow(index.mat), ncol=ncol(index.mat), dimnames=dimnames(index.mat))
    
    corrplot(comparison.mat, method="square", is.corr=FALSE, cl.lim=c(-1,1), tl.srt=45, col=cols(200), tl.col="black")
    corrplot(comparison.mat, method="square", is.corr=FALSE, cl.lim=c(-1,1), tl.srt=45, bg=comparison.cuts, add=TRUE, addgrid.col=NA, col=alpha("white",0), cl.pos="n", tl.pos="n")
    corrplot(index.mat, method="square", is.corr=FALSE, cl.lim=c(-1,1), tl.srt=45, bg=alpha("white",0), col=index.cuts, add=TRUE, addgrid.col=NA, cl.pos="n", , tl.pos="n", p.mat=index.mat.p, sig.level=.05)
    text(-1.3, 9.5, "Legend", font=2)
    rect(-1, 8.5, 0, 9.5, col=cols(200)[150], border=NA)
    rect(-.8, 8.7, -.2, 9.3, col=cols(200)[120], border=NA)
    text(-.5, 9.4, "present")
    text(-.5, 9, "track\nindex")
}

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
pdf("figures/traits-05.pdf"); .plot.corr.mat(stem, "b.track.index", "pres.dis.pres.env", "0.05", "observed", clean=5); dev.off()
pdf("figures/traits-25.pdf"); .plot.corr.mat(stem, "b.track.index", "pres.dis.pres.env", "0.25", "observed", 5); dev.off()
pdf("figures/traits-50.pdf"); .plot.corr.mat(stem, "b.track.index", "pres.dis.pres.env", "0.5", "observed", 5); dev.off()
pdf("figures/traits-75.pdf"); .plot.corr.mat(stem, "b.track.index", "pres.dis.pres.env", "0.75", "observed", 5); dev.off()
pdf("figures/traits-95.pdf"); .plot.corr.mat(stem, "b.track.index", "pres.dis.pres.env", "0.95", "observed", 5); dev.off()

# Correlations
.get.vals <- function(stem, index, null, clean, p)
    return(
        abind(
            .corr.mat(stem, index, "0.05", null, clean, p),
            .corr.mat(stem, index, "0.25", null, clean, p),
            .corr.mat(stem, index, "0.5", null, clean, p),
            .corr.mat(stem, index, "0.75", null, clean, p),
            .corr.mat(stem, index, "0.95", null, clean, p),
            along=3
        )
    )
index.r <- .get.vals(stem, "b.track.index", "observed", 100, FALSE)
index.p <- .get.vals(stem, "b.track.index", "observed", 100, TRUE)
past.r <- .get.vals(stem, "past.dis.past.env", "observed", 100, FALSE)
past.p <- .get.vals(stem, "past.dis.past.env", "observed", 100, TRUE)
pres.r <- .get.vals(stem, "pres.dis.pres.env", "observed", 100, FALSE)
pres.p <- .get.vals(stem, "pres.dis.pres.env", "observed", 100, TRUE)
sum(index.p < .05); prod(dim(index.p))
sum(past.p < .05); prod(dim(index.p))
sum(pres.p < .05); prod(dim(index.p))



plot(as.numeric(index.r) ~ as.numeric(past.r), asp=1, xlab="", ylab="")
abline(0, 1, col="grey40", lty=2, lwd=3)


# Making the combined traits and phylogeny plot
# Get data
traits <- data.frame(index=as.numeric(index.r[,,3]), past=as.numeric(past.r[,,3]), pres=as.numeric(pres.r[,,3]), trait=rep(rownames(index.r),ncol(index.r)), climate=rep(colnames(index.r),each=nrow(index.r)), taxon=rep(c("mammals","birds","plants","plants","plants","plants"),ncol(index.r)), type="traits")
phylo <- data.frame(
    index=c(
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/birds-signal.RDS")$bootstrap.index["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/mammals-signal.RDS")$bootstrap.index["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/plants-signal.RDS")$bootstrap.index["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/reptiles-signal.RDS")$bootstrap.index["0.5",])
    ),
    past=c(
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/birds-signal.RDS")$past["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/mammals-signal.RDS")$past["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/plants-signal.RDS")$past["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/reptiles-signal.RDS")$past["0.5",])
    ),
    pres=c(
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/birds-signal.RDS")$present["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/mammals-signal.RDS")$present["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/plants-signal.RDS")$present["0.5",]),
        as.numeric(readRDS("../track-index-2019-r1+/clean-data/reptiles-signal.RDS")$present["0.5",])
    ),
    trait="phylo",
    climate=rep(colnames(readRDS("../track-index-2019-r1+/clean-data/reptiles-signal.RDS")$bootstrap.index), 4),
    taxon=rep(c("mammals","birds","plants","reptiles"), each=9),
    type="phylo"
)
data <- rbind(traits, phylo)

z.off <- .1
data$index <- abs(data$index)+z.off; data$pres <- abs(data$pres)+z.off; data$past <- abs(data$past)+z.off
data$index[data$type=="phylo"] <- -data$index[data$type=="phylo"]; data$pres[data$type=="phylo"] <- -data$pres[data$type=="phylo"]; data$past[data$type=="phylo"] <- -data$past[data$type=="phylo"]
data$t.taxon <- as.numeric(factor(data$taxon))
data$t.trait <- as.numeric(factor(data$trait))
lookup <- c("clouds"="cld", "frost"="frs", "evapo-trans."="vap", "precipitation"="pre", "min(temp)"="tmn", "mean(temp)"="tmn", "max(temp)"="tmx", "vapour"="vap", "rainy day"="wet", "cld"="cld", "frs"="frs", "vap"="vap", "pre"="pre", "tmn"="tmn", "tmp"="tmp", "tmx"="tmx", "vap"="vap", "wet"="wet", "pet"="pet")
data$climate <- lookup[data$climate]
data$t.climate <- (as.numeric(factor(data$climate))-4.5)/(4.5*4)

pdf("figures/trait-phylo.pdf")
with(data, plot(index ~ I(t.taxon+t.climate), pch=ifelse(abs(index)>abs(past) & abs(index)>abs(pres),20,21), xlab="", ylab="", axes=FALSE, ylim=c(-1-z.off,1+z.off), col="black"))
with(data, points(pres ~ I(t.taxon+t.climate), pch=ifelse(abs(pres)>abs(past) & abs(pres)>abs(index),20,21), col="blue"))
with(data, points(past ~ I(t.taxon+t.climate), pch=ifelse(abs(past)>abs(index) & abs(past)>abs(pres),20,21), col="red"))
axis(2, at=seq(z.off,1+z.off,by=.2),  labels=c(NA, seq(.2,1,by=.2)))
axis(2, at=-seq(z.off,1+z.off,by=.2), labels=c(NA, seq(.2,1,by=.2)))
#abline(h=0, col="grey40", lwd=5)
mtext("trait correlations", side=2, adj=.5, at=.5, line=2, cex=1.2)
mtext("phylogenetic signal", side=2, adj=.5, at=-.5, line=2, cex=1.2)
legend(.8, 1.2, pch=20, col=c("black","blue","red"), legend=c("index","past","present"), cex=1.2, bty="n")
legend(1.6, 1.2, pch=c(20,21), col="black", legend=c("strongest","weaker"), cex=1.2, bty="n")
grid.raster(readPNG("phylopics/mammal.png"),   .20, .52, width=.06)
grid.raster(readPNG("phylopics/bird.png"),     .42, .52, width=.06)
grid.raster(readPNG("phylopics/plant.png"),    .65, .52, width=.06)
grid.raster(readPNG("phylopics/reptile.png"),  .88, .50, width=.06)
dev.off()

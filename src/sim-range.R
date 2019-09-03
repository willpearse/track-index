# Load data
source("src/headers.R")
temp <- readRDS("clean-data/cru-tmp.RDS")

# Wrapper functions
.make.block.dist <- function(x, y, range)
    as.matrix(expand.grid(x+seq(-range,range), y+seq(-range,range)))
.check.viable <- function(dist, raster)
    !(any(dist < 0) | any(dist[,1] > nrow(raster)) | any(dist[,2] > ncol(raster)))
.get.niche <- function(dist, raster)
    median(raster[dist], na.rm=TRUE)
.prob.occur <- function(dist, raster, alpha, sd, optimum){
    env <- raster[dist]
    niche <- dnorm(env, optimum, sd) / dnorm(optimum, optimum, sd)
    niche <- (1-alpha) + niche*alpha
    present <- dist[runif(length(env)) <= niche,]
    return(present)
}
.dist.wrap <- function(x, y, range, raster, alpha, sd, optimum){
    block <- .make.block.dist(x, y, range)
    return(.prob.occur(block, raster, alpha, sd, optimum))
}
.do.calc <- function(i){
    # Make distributions
    optimum <- .get.niche(.make.block.dist(data$x[i], data$y[i], data$range[i]), early)
    early.dist <- .dist.wrap(data$x[i], data$y[i], data$range[i], early, data$alpha[i], 1, optimum)
    mid.dist <- .dist.wrap(data$x[i], data$y[i], data$range[i], mid, data$alpha[i], 1, optimum)
    late.dist <- .dist.wrap(data$x[i], data$y[i], data$range[i], late, data$alpha[i], 1, optimum)
    # Calculate index
    mid.index <- .calc.metric(as.numeric(mid[mid.dist]), as.numeric(early[early.dist]), as.numeric(mid[early.dist]), quantile=.5, n.boot=99)
    late.index <- .calc.metric(as.numeric(late[late.dist]), as.numeric(early[early.dist]), as.numeric(late[early.dist]), quantile=.5, n.boot=99)
    # Check for successful completion and return
    if(is.null(mid.index) | is.null(late.index))
        return(matrix(NA, 2, 7, dimnames=list(c("mid","late"), c("index","bootstrap.index", "mad.index", "quantile","present","past","projected"))))
    return(matrix(c(mid.index,late.index), 2, 7, dimnames=list(c("mid","late"), c("index","bootstrap.index", "mad.index", "quantile","present","past","projected")), byrow=TRUE))
}

# Prepare pre- and post- data
early <- as.matrix(temp[[1]])
mid <- as.matrix(temp[[50]])
late <- as.matrix(temp[[100]])

# Select viable starting distributions
centers <- which(!is.na(as.matrix(early)), arr.ind=TRUE)
centers <- centers[sample(nrow(centers), 100),]
viable <- rep(FALSE, nrow(centers))
for(i in seq_along(viable))
    viable[i] <- .check.viable(.make.block.dist(centers[i,1], centers[i,2], 10), early)
centers <- centers[viable,]

# Add parameters for simulations
species <- expand.grid(shift=seq(-2,2), range=c(1,2,5), alpha=0, merge=paste(centers[,1], centers[,1], sep="-"))
species <- rbind(species,
                 expand.grid(shift=seq(-2,2), range=c(5,10,20,50), alpha=c(0,.5,1), merge=paste(centers[,1], centers[,1], sep="-")))
centers <- data.frame(x=centers[,1], y=centers[,2], merge=paste(centers[,1], centers[,1], sep="-"))
data <- merge(species, centers, by="merge")

# Do work
results <- mcMap(.do.calc, seq_len(nrow(data)))
save.image("sim-range.RData")


# Some post-processing
for(i in seq_along(late.index)){
    if(is.null(late.index[[i]]))
        late.index[[i]] <- setNames(rep(NA, 7), c("index","bootstrap.index", "mad.index", "quantile","present","past","projected"))
    if(is.null(mid.index[[i]]))
        mid.index[[i]] <- setNames(rep(NA, 7), c("index","bootstrap.index", "mad.index", "quantile","present","past","projected"))
}
late.index <- do.call(rbind, late.index)
late.index <- cbind(data, late.index)
mid.index <- do.call(rbind, mid.index)
mid.index <- cbind(data, mid.index)
late.index <- na.omit(late.index)
mid.index <- na.omit(mid.index)

# Assign what happened
mid.index$outcome <- with(mid.index, mapply(.assign.outcome, present, past, projected))
late.index$outcome <- with(late.index, mapply(.assign.outcome, present, past, projected))
outcome.col <- setNames(c("forestgreen","blue","grey60","red"), c("lag","overshoot","static-distribution","undershoot"))

# "Remove" outliers
mid.index$backup <- mid.index$index
late.index$backup <- late.index$index
mid.index$index[abs(mid.index$index) > 50] <- NA
late.index$index[abs(late.index$index) > 50] <- NA

#mid.index$jitter <- jitter(mid.index$shift)
mid.index$c.change <- with(mid.index, ifelse(present-past > 0, "red", "blue"))
mid.index$c.change[mid.index$shift==0] <- "black"

mid.index <- mid.index[order(mid.index$shift, mid.index$index),]
unique.shifts <- sort(unique(mid.index$shift))
mid.index$jitter <- -1
for(i in seq_along(unique.shifts)){
    which <- which(mid.index$shift==unique.shifts[i])
    mid.index$jitter[which] <- mid.index$shift[which] + seq(-.4,.4, length.out=length(which))
}

pdf("null-movement.pdf", width=10, height=5)
with(mid.index, plot(bootstrap.index ~ jitter, ylim=c(-10,10), pch=20, col=outcome.col[outcome], axes=FALSE, xlab="Grid cell shift", ylab="Bootstrapped index"))
axis(1)
axis(2)
with(mid.index, arrows(x0=jitter, y0=bootstrap.index-mad.index, y1=bootstrap.index+mad.index, col=outcome.col[outcome], lwd=.5, code=3, angle=0))
abline(h=1, col="grey60", lwd=3)
abline(h=0, col="forestgreen", lwd=3)
legend("topright", pch=20, col=outcome.col, legend=names(outcome.col), bty="n")
dev.off()

with(mid.index, summary(lm(index ~ outcome -1)))

#######################################################################
## readcel.R for spatial correction
## 2007-09-06 copied from natural.uchicago.edu's /n1/naturalsystems/array/part2/atsnptile1/readcel.R
library(affy)
library(stats)


readcel <- function(cel.files, probe.number, array.size, filter.size ){ 

mprobe.mean <- matrix(NA, nr=probe.number, nc=length(cel.files) )
for (i in 1:length(cel.files) )
{
probe.mean <- log(intensity(read.affybatch (filenames=cel.files [i] ) ) )
mprobe.mean[,i] <- matrix(probe.mean, nc=array.size, byrow=T) [array.size:1,] [probe.ok] 
}

# calculate median probe within array
meds <- apply(mprobe.mean,2,median)
mprobe.resid <- mprobe.mean - meds
filter.size <- 2*trunc(filter.size /2)+1  # force to be odd so the edge will be same length = filter.size/2
fe <- (filter.size-1)/2  # filter edge
mprobe.bg <- mprobe.resid

# do spatial correction
for (i in 1:length(cel.files))
{
    sp.mat <- matrix(0,nr= array.size,nc=array.size)
    sp.mat[probe.ok] <- mprobe.resid[,i]
    sp.mat.sum <- apply(sp.mat,2,filter,rep(1,filter.size))  # from start to end, every 51 probes sum in moving 
    sp.mat.sum <- apply(sp.mat.sum,1,filter,rep(1,filter.size))
    sp.mat.count <- apply(probe.ok,2,filter,rep(1,filter.size)) 
    sp.mat.count <- apply(sp.mat.count,1,filter,rep(1,filter.size))  #how many probes have been added up to contribute the sp.mat.sum
    spatial.corrected <- t(sp.mat.sum/sp.mat.count)

    # fill in edges, with the background of edge boarder
    spatial.corrected[1:fe,] <- rep(spatial.corrected[fe+1, ], each = fe)
    spatial.corrected[(array.size+1 - fe):array.size, ] <-
                        rep(spatial.corrected[array.size - fe, ], each = fe)
    spatial.corrected[,1:(fe)] <- spatial.corrected[,fe+1]
    spatial.corrected[,(array.size+1 - fe):array.size] <-
                                      spatial.corrected[,array.size - fe]
    
    mprobe.bg[,i] <- spatial.corrected[probe.ok]
}

 mprobe.mean <- (mprobe.mean - mprobe.bg)[chrom.order,]
 return (mprobe.mean)

}

#######################################################################
## function for reading in cel files and performing spatial correction
## 2007-09-06 copied from natural.uchicago.edu's /n1/naturalsystems/array/part2/atsnptile1/readcel2.R
readcel2 <- function(cel.files, array.size = 712, filter.size = 51, imageType="jpg", probe.OK=""){
# read .cel files

if (array.size > 720){ # for tiling arrays
filter.size <- 81
mprobe.mean <- matrix(ncol = length(cel.files),nrow = length(chrom.order))
for (i in 1:length(cel.files))
   ## turn image in by rotate right 90 so xy coords match probe.OK
   mprobe.mean[,i] <- log( matrix( intensity( read.affybatch(filenames = cel.files[i])),
              nc = array.size, byrow = T)[array.size:1,][probe.OK])
}

#if (array.size != 2560){ 
#mprobe.mean <- log(intensity(ReadAffy(filenames=cel.files))[c(probe.OK),])
#}

# calculate median probe across arrays
meds <- apply(mprobe.mean,1,median)
mprobe.resid <- mprobe.mean - meds
filter.size <- 2*trunc(filter.size /2)+1  ## force to be odd
fe <- (filter.size-1)/2 # filter edge
mprobe.bg <- mprobe.resid
# do spatial correction
if(imageType == "pdf") pdf(file = "spatialCorrected.pdf")
for (i in 1:length(cel.files)){
    sp.mat <- matrix(0,nr= array.size,nc=array.size)
    sp.mat[probe.OK] <- mprobe.resid[,i]
    sp.mat.sum <- apply(sp.mat,2,filter,rep(1,filter.size))	#column-wise filter
    sp.mat.sum <- apply(sp.mat.sum,1,filter,rep(1,filter.size))	#row-wise filter
    sp.mat.count <- apply(probe.OK,2,filter,rep(1,filter.size))	#column-wise 
    sp.mat.count <- apply(sp.mat.count,1,filter,rep(1,filter.size))
    spatial.corrected <- t(sp.mat.sum/sp.mat.count)	#?
    ## fill in edges
    spatial.corrected[1:fe,] <- rep(spatial.corrected[fe+1, ], each = fe)
    spatial.corrected[(array.size+1 - fe):array.size, ] <-
                        rep(spatial.corrected[array.size - fe, ], each = fe)
    spatial.corrected[,1:(fe)] <- spatial.corrected[,fe+1]
    spatial.corrected[,(array.size+1 - fe):array.size] <-
                                      spatial.corrected[,array.size - fe]
    sp.range <- round(range(spatial.corrected),1)
    if (imageType=="jpg"){
        jpeg(file=paste(cel.files[i],"sc.jpg",sep=""))
        image(1:array.size, 1:array.size, spatial.corrected, col = terrain.colors(2^8),
            main = paste(cel.files[i], "spatialResid", filter.size,max(spatial.corrected),min(spatial.corrected),sep=".") )
        dev.off()
    }
    if (imageType=="pdf") image(1:array.size, 1:array.size, spatial.corrected, col = terrain.colors(2^8),
                  main = paste(cel.files[i], "spatialResid", filter.size,sp.range,sep="") )
    mprobe.bg[,i] <- spatial.corrected[probe.OK]
    }
if(imageType == "pdf") dev.off()
(mprobe.mean - mprobe.bg)[chrom.order,]
}

## 2007-09-06 copied form Xu Zhang's email
load ("snp.RData")
library(affy)
cel.files <- list.celfiles()


ord <- rep( c(1,2,3,4), nrow(snp)/4 )
snp <- cbind( snp, ord)
array.size <- 1612
probe.ok <- matrix(FALSE, nr=array.size, nc=array.size)
probe.ok.chrom <- matrix(NA, nr=array.size, nc=array.size)
probe.ok.position <- matrix(NA, nr=array.size, nc=array.size)
probe.ok.ord <- matrix(NA, nr=array.size, nc=array.size)
probe.ok[cbind(snp$xpos+1, snp$ypos+1)] <- TRUE
probe.ok.chrom[cbind(snp$xpos + 1, snp$ypos+1)] <-snp$chr
probe.ok.position [cbind(snp$xpos + 1, snp$ypos+1)] <- snp$bp
probe.ok.ord [cbind(snp$xpos + 1, snp$ypos+1)] <- snp$ord

probe.ok.chrom.valid <- probe.ok.chrom[probe.ok]
probe.ok.position.valid <- probe.ok.position [probe.ok]
probe.ok.ord.valid <- probe.ok.ord [probe.ok]
chrom.order <- order(probe.ok.chrom.valid, probe.ok.position.valid, probe.ok.ord.valid)


mprobe.mean <- readcel (cel.files=cel.files, probe.number=nrow(snp), array.size=1612, filter.size=81)

mprobe.mean <- normalize.quantiles(mprobe.mean)
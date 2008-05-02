### readcel.R for spatial correction
library(affy)
library(stats)
readcel <- function(cel.files, probe.number, array.size, filter.size ){ 

mprobe.mean <- matrix(NA, nr=probe.number, nc=length(cel.files) )
for (i in 1:length(cel.files) )
{
probe.mean <- log(intensity(read.affybatch(filenames=cel.files [i] ) ) )
#modified by R.J. below
#probe.mean <- log(intensity(ReadAffy(filenames=cel.files [i] ) ) )


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
    sp.mat[!probe.ok] <- median(mprobe.resid[,i])

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

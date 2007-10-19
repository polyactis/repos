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
readcel2 <- function(cel_file_directory, array.size = 712, filter.size = 51, imageType="jpg", probe.OK=""){
# read .cel files
cel.files <- list.celfiles(cel_file_directory, full.names=TRUE)

if (array.size > 720){ # for tiling arrays
filter.size <- 81
mprobe.mean <- matrix(ncol = length(cel.files),nrow = length(chrom.order))
for (i in 1:length(cel.files))
{
   cat(cel.files[i])
   ## turn image in by rotate right 90 so xy coords match probe.OK
   mprobe.mean[,i] <- log( matrix( intensity( read.affybatch(filenames = cel.files[i])),
              nc = array.size, byrow = T)[array.size:1,][probe.OK])
   cat(".\n")
}
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
cat("spatial correction\n")
for (i in 1:length(cel.files)){
    cat("\t", cel.files[i])
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
    cat(".\n")
}

if(imageType == "pdf") dev.off()

#(mprobe.mean - mprobe.bg)[chrom.order,]
mprobe.mean[chrom.order,]
}

#######################################################################
##2007-09-06 function to read raw samples from a directory to make a matrix, no log transformation
readcel_raw <- function(cel_file_directory, array.size = 712, probe.OK="", chrom.order=""){
# read .cel files
cel.files <- list.celfiles(cel_file_directory, full.names=TRUE)

mprobe.mean <- matrix(ncol = length(cel.files),nrow = length(chrom.order))
for (i in 1:length(cel.files))
{
   cat(cel.files[i])
   ## turn image in by rotate right 90 so xy coords match probe.OK
   mprobe.mean[,i] <- matrix( intensity( read.affybatch(filenames = cel.files[i])),
              nc = array.size, byrow = T)[array.size:1,][probe.OK]
   cat(".\n")
}
mprobe.mean[chrom.order,]
}


#######################################################################
##2007-09-07 function to check sign of mprobe.raw data, 1 means pairs of two sense strains and two anti-sense are in same order (either > or <). -1 means in discordant order. 0 means one of the pairs show equal intensity count.
check_sign <- function(data)
{
no_of_snps=dim(data)[1]/4
no_of_arrays=dim(data)[2]
for (j in 1:no_of_arrays)
{
  sign_ls = matrix(0, nrow=1, ncol=no_of_snps)
  d1_ls = matrix(0, nrow=1, ncol=no_of_snps)
  d2_ls = matrix(0, nrow=1, ncol=no_of_snps)
  for (i in 1:no_of_snps)
  {
    i = i-1	#then i starts from 0
    d1 = data[i*4+1,j]-data[i*4+3,j]
    d2 = data[i*4+2,j]-data[i*4+4,j]
    d3 = log(data[i*4+1,j]/data[i*4+2,j],2)
    d4 = log(data[i*4+3,j]/data[i*4+4,j],2)
    d5 = 1/4*(log(data[i*4+1,j])+log(data[i*4+2,j])+log(data[i*4+3,j])+log(data[i*4+4,j]))
    d1_ls[i+1] = (d3+d4)/2
    d2_ls[i+1] = d5
    sign_ls[i+1] = sign(d1*d2)
  }
  sign_count=table(sign_ls)
  print(sign_count)
  print(sign_count[1]/no_of_snps)
  plot(d2_ls, d1_ls)
}
}

#######################################################################
##2007-09-07 function to draw a MA-like plot. it's only for 6 arrays from 'yanli_8-8-07'. function readcel_raw() is used to read these 6 arrays. each strain got two copies. so a simple MA-plot which compares two arrays' intensity
ma_p <- function(data)
{
no_of_snps=dim(data)[1]
no_of_arrays=dim(data)[2]
for (j in 1:as.integer(no_of_arrays/2))
{
  sign_ls = matrix(0, nrow=1, ncol=no_of_snps)
  d1_ls = matrix(0, nrow=1, ncol=no_of_snps)
  d2_ls = matrix(0, nrow=1, ncol=no_of_snps)
  j = j-1
  for (i in 1:no_of_snps)
  {
    d1 = log2(data[i,j*2+1]/data[i,j*2+2])
    d2 = 1/2*log2(data[i,j*2+1]*data[i,j*2+2])
    d1_ls[i] = d1
    d2_ls[i] = d2
    sign_ls[i] = sign(d1*d2)
  }
  sign_count=table(sign_ls)
  print(sign_count)
  print(sign_count[1]/no_of_snps)
  plot(d2_ls, d1_ls)
}
}


## 2007-09-06 code starts from here. based on snippets copied form Xu Zhang's email
load ("snp.RData")
library(affy)

##data frame 'snp' is loaded from 'snp.RData', which has all sorts of information about snp probes. check 'snp' is already loaded in by objects().
objects()
##print the 1st 10 probes
snp[1:10,]
ord <- rep( c(1,2,3,4), nrow(snp)/4 )	##every row in 'snp' is a probe for one snp. there're 4 probes for one snp, order them 1,2,3,4.
snp <- cbind( snp, ord)

##2007-09-07 convert some columns into integer class to avoid the bug. this is the way i came up. probably there's a better way to do it.
write.table(snp, '/tmp/snp', sep='\t')
snp=read.table('/tmp/snp',as.is=c(1:4,7:9))	#column 1,2,3,4,7,8,9 is snpid, seq, chr, bp, xpos, ypos, ord

##another way to do the conversoin: snp[,4] = as.integer(as.character(snp[,4]))

##map most useful columns of 'snp' to a matrix whose cells correspond to that of the physical array
array.size <- 1612
probe.ok <- matrix(FALSE, nr=array.size, nc=array.size)
probe.ok.chrom <- matrix(NA, nr=array.size, nc=array.size)
probe.ok.position <- matrix(NA, nr=array.size, nc=array.size)
probe.ok.ord <- matrix(NA, nr=array.size, nc=array.size)
probe.ok[cbind(snp$xpos+1, snp$ypos+1)] <- TRUE
probe.ok.chrom[cbind(snp$xpos + 1, snp$ypos+1)] <-snp$chr
probe.ok.position [cbind(snp$xpos + 1, snp$ypos+1)] <- snp$bp	#when factor is converted into integer, it messed up. 657 => 202179
probe.ok.ord [cbind(snp$xpos + 1, snp$ypos+1)] <- snp$ord
probe.ok.chrom.valid <- probe.ok.chrom[probe.ok]
probe.ok.position.valid <- probe.ok.position [probe.ok]
probe.ok.ord.valid <- probe.ok.ord [probe.ok]
chrom.order <- order(probe.ok.chrom.valid, probe.ok.position.valid, probe.ok.ord.valid)

##show how probes look like in chrom.order
snp_chr_order_data = rbind(probe.ok.chrom.valid[chrom.order], probe.ok.position.valid[chrom.order], probe.ok.ord.valid[chrom.order])
snp_chr_order_data[,1:20]

##now start to read files, readcel2() did log transformation, filtering, etc.
cel_file_directory = 'yanli_8-8-07'
mprobe.mean <- readcel2('yanli_8-8-07', array.size=1612, filter.size=81, imageType="jpg", probe.OK=probe.ok)
mprobe.mean <- normalize.quantiles(mprobe.mean)

##read the raw intensity data
mprobe.raw = readcel_raw(cel_file_directory, array.size = 1612, probe.OK=probe.ok, chrom.order=chrom.order)

##2007-09-09 show where each chromsome's probes are
xls=1:array.size
#chromsome 3
image(xls, 1:array.size, (probe.ok.chrom==3), col=topo.colors(2^8))
#all chromsomes
jpeg(file='genotyping_probes_on_the_array.jpeg')
image(xls, 1:array.size, (probe.ok.chrom==3)+(probe.ok.chrom==4)+(probe.ok.chrom==1)+(probe.ok.chrom==2)+(probe.ok.chrom==5), col=topo.colors(2^8), main='location of genotyping probes on the array')
dev.off()
##a simple sign test to see how 
check_sign(mprobe.raw)

##MA-like plot for 6 replicate arrays in 'yanli_8-8-07'
ma_p(mprobe.raw)

##2007-09-07 just test, on the way to find out the bug
for (i in 1:dim(probe.ok.position)[1])
   for (j in 1:dim(probe.ok.position)[2])
     if (!is.na(probe.ok.position[i,j]))
     {
	if(probe.ok.position[i,j]==202179)
        {
          print(probe.ok.chrom[i,j])
          print(probe.ok.ord[i,j])
        }
      }
### filename: new_basecall.R
### content: make base calls from raw cel file, namely "raw_data.cel"
###  output: genocall.txt where 1 is the first allele and 3 is the second allele
###          the 2nd column is the likelihood of correct calls
###    note: use  to transform genocall.txt to tsv file given a cutoff on the likelihood
###	
###	R. J. 04/25/2008

source("utilities-jr.R")
source ("readcel.R")
source("cal-pacc.R")

library(splines);
library(affy)

load ("snp.RData")  #annotation for SNP probes
load("params-logitreg-31arrays.rda")
load("params-31arrays.rda")
load("m.norm.107arrays.rda");
m.norm <- m.norm/107;


#added specification for the Single CEL file
cel.files <- "raw_data.cel";



#read in cel intensities
ord <- rep( c(1,2,3,4), nrow(snp)/4 )
snp <- cbind( snp, ord)
array.size <- 1612
probe.ok <- matrix(FALSE, nr=array.size, nc=array.size)
probe.ok.chrom <- matrix(NA, nr=array.size, nc=array.size)
probe.ok.position <- matrix(NA, nr=array.size, nc=array.size)
probe.ok.ord <- matrix(NA, nr=array.size, nc=array.size)
probe.ok[cbind(snp$xpos+1, snp$ypos+1)] <- TRUE
probe.ok.chrom[cbind(snp$xpos + 1, snp$ypos+1)] <- as.numeric( as.character(snp$chr) )
probe.ok.position [cbind(snp$xpos + 1, snp$ypos+1)] <- as.numeric( as.character( snp$bp) )
probe.ok.ord [cbind(snp$xpos + 1, snp$ypos+1)] <- as.numeric( as.character(snp$ord) )
 
probe.ok.chrom.valid <- probe.ok.chrom[probe.ok]
probe.ok.position.valid <- probe.ok.position [probe.ok]
probe.ok.ord.valid <- probe.ok.ord [probe.ok]
chrom.order <- order(probe.ok.chrom.valid, probe.ok.position.valid, probe.ok.ord.valid)


mprobe.mean <- readcel (cel.files=cel.files, probe.number=nrow(snp), array.size=1612, filter.size=81) #background substraction

#quantile normalization
#m.norm is the average among sorted mprobe.mean of 107 arrays
tmprank <- rank(mprobe.mean);
new.m.norm <- (m.norm*107 + mprobe.mean)/108
mprobe.norm <- new.m.norm[tmprank];
save(mprobe.norm, file="mprobe_norm.rda");

#call genotypes
dat <- matrix(mprobe.norm, ncol = 1);
calls <- genocall(dat);
pacc <- calpacc(dat, index = 1:(nrow(dat)/4), calls$ff, calls$LLR);
tmp <- cbind(calls$myCalls, pacc);
write.table(tmp, file= "genocalls.txt", quote=FALSE, row.names=FALSE, col.names=FALSE);

q("no");

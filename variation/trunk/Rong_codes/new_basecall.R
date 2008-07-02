#!/usr/bin/env Rscript
### filename: new_basecall.R
### content: make base calls from raw cel file, namely "raw_data.cel"
###  output: genocall.txt where 1 is the first allele and 3 is the second allele
###          the 2nd column is the likelihood of correct calls
###    note: use  to transform genocall.txt to tsv file given a cutoff on the likelihood
###	
###	R. J. 04/25/2008
command_args = commandArgs()
program_name = substr(command_args[4], 8, nchar(command_args[4]))	#command_args[4] = "--file=./new_basecall.R"

if (length(command_args)!=9)	#2008-07-01 command_args starts with ["/usr/lib/R/bin/exec/R", "--slave", "--no-restore", "--file=../test/R/test_command_arguments.R", "--args"]
{
	cat("\n")
	cat("Usage:", program_name, "INPUT_DIR START_ARRAY_ID END_ARRAY_ID OUTPUT_DIR\n")
	cat("\n")
	cat("	Do oligo genotype call on arrays in input_dir(/Network/Data/250k/db/raw_data/) \n")
	cat("	it loads some other R scripts and R data file in the same directory. It assumes cel files are named like $(array_id)_raw_data.cel in the input directory.\n")
	q('no')	#dont' save this session
}

input_dir=command_args[6]
start_array_id=as.integer(command_args[7])
end_array_id=as.integer(command_args[8])
output_dir=command_args[9]

cat("Loading data ...")
source("utilities-jr.R")
source("readcel.R")
source("cal-pacc.R")

library(splines);
library(affy)

load ("snp.RData")  #annotation for SNP probes
load("params-logitreg-31arrays.rda")
load("params-31arrays.rda")
load("m.norm.107arrays.rda");
m.norm <- m.norm/107;
cat("\n")

cat("Preparing probe data structure ...")
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
cat("\n")

for (i in seq(start_array_id, end_array_id))
{
	#added specification for the Single CEL file
	cel.files <- paste(input_dir, '/', i , '_raw_data.cel', sep="")
	cat("Array: ", i, cel.files, "...")
	
	mprobe.mean <- readcel (cel.files=cel.files, probe.number=nrow(snp), array.size=1612, filter.size=81) #background substraction
	
	#quantile normalization
	#m.norm is the average among sorted mprobe.mean of 107 arrays
	tmprank <- rank(mprobe.mean);
	new.m.norm <- (m.norm*107 + mprobe.mean)/108
	mprobe.norm <- new.m.norm[tmprank];
	save(mprobe.mean, file=paste(output_dir, '/', i , '_mprobe_mean.rda', sep=""));
	save(mprobe.norm, file=paste(output_dir, '/', i , '_mprobe_norm.rda', sep=""));
	
	#call genotypes
	dat <- matrix(mprobe.norm, ncol = 1);
	calls <- genocall(dat);
	pacc <- calpacc(dat, index = 1:(nrow(dat)/4), calls$ff, calls$LLR);
	tmp <- cbind(calls$myCalls, pacc);
	write.table(tmp, file=paste(output_dir, '/', i , '_genocalls.txt', sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE);
	cat("\n")
}

q("no");	#dont' save this session

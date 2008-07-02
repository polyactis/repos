#!/usr/bin/env Rscript
#2008-05-17 Rong's piece of code to draw histogram of normalized intensity and unnormalized intensity from each array.
#modifed to take arguments, for loop, etc
#R --vanilla --args $0 $* << EOF	#old way to make R script standalone

command_args = commandArgs()
#program_name = command_args[4]
program_name = substr(command_args[4], 8, nchar(command_args[4]))	#command_args[4] = "--file=./plot-mean-...R"

if (length(command_args)!=9)	#2008-07-01 command_args starts with ["/usr/lib/R/bin/exec/R", "--slave", "--no-restore", "--file=../test/R/test_command_arguments.R", "--args"]
{
	cat("\n")
	cat("Usage:", program_name, "INPUT_DIR START_ARRAY_ID END_ARRAY_ID PDF_OUTPUT_FNAME\n")
	cat("\n")
	cat("	Plot histogram of normalized intensity and unnormalized intensity from each array\n")
	q('no')	#dont' save this session
}

input_dir = command_args[6]
start_array_id = as.real(command_args[7])
end_array_id = as.real(command_args[8])
output_fname = command_args[9]

load("m.norm.107arrays.rda")	#take this as normalized data
m.norm <- m.norm/107
pdf(output_fname)

for (i in seq(start_array_id, end_array_id))
{
	cat("Array:", i)
	#1. plot the learning set, 107 arrays' intensity
	plot(density(m.norm), "l", xlim=c(1, 10), ylim=c(0, 1), main=paste("array id = ", i))
	
	#2. plot the original intensity of this array
	probe_mean_fname = paste(input_dir, '/', i , '_mprobe_mean.rda', sep="")	
	load(probe_mean_fname)
	lines(density(mprobe.mean), col="red")
	
	#3. plot the intensity after normalization
	probe_norm_fname = paste(input_dir, '/', i , '_mprobe_norm.rda', sep="")
	load(probe_norm_fname)
	lines(density(mprobe.norm), col="blue")
	
	rm(mprobe.norm)
	rm(mprobe.mean)
	cat("\n")
}

dev.off()
q("no")

EOF

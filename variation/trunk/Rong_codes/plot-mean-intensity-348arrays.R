#2008-05-17 Rong's piece of code to draw histogram of normalized intensity and unnormalized intensity from each array.
#modifed to take arguments, for loop, etc
R --vanilla --args $0 $* << EOF

command_args = commandArgs()
print(command_args)
program_name = command_args[4]

if (length(command_args)!=8)	#2008-05-17 command_args starts with ["/usr/lib/R/bin/exec/R", "--vanilla", "--args"]
{
	cat("\n")
	cat("Usage:", program_name, "INPUT_DIR START_ARRAY_ID END_ARRAY_ID PDF_OUTPUT_FNAME\n")
	cat("\n")
	cat("	Plot histogram of normalized intensity and unnormalized intensity from each array\n")
	q('no')	#dont' save this session
}

input_dir = command_args[5]
start_array_id = as.real(command_args[6])
end_array_id = as.real(command_args[7])
output_fname = command_args[8]

load("m.norm.107arrays.rda")	#take this as normalized data
m.norm <- m.norm/107
pdf(output_fname)

for (i in seq(start_array_id, end_array_id))
{
	probe_mean_fname = paste(input_dir, '/', i , '_mprobe_mean.rda', sep="")
	cat(probe_mean_fname)
	load(probe_mean_fname)
	plot(density(m.norm), "l", xlim=c(1, 10), ylim=c(0, 1), main=paste("array id = ", i))
	lines(density(mprobe.mean), col="red")
	rm(mprobe.mean)
	cat("\n")
}

dev.off()
q("no")

EOF

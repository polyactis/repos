%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<BEAM 1.0>>
Yu Zhang, Jun S Liu

This is a brief guide for using BEAM (Bayesian Epistatis Association Mapping). 
The program uses Markov Chain Monte Carlo (MCMC) to search for both single-marker
and interaction effects from case-control SNP data.

Reference:

Zhang Y and Liu JS (2007). Bayesian Inference of Epistatic Interations in Case-Control Studies.
Nature Genetics, in press.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This package contains the following files:	
	BEAM
	libgsl.dll  		<-- for dos version only
	libgslcblas.dll		<-- for dos version only
	parameters.txt
	README.txt

[1. Installation]-----------------------------------------

Unzip all files into one folder. BEAM uses the GNU Scientific Library (GSL), for 
which two necessary files "libgsl.dll" and "libslcblas.dll" are included in
the DOS package. Make sure you put the two dll files in the same folder with BEAM,
or put them in your system folder (e.g., "c:\\windows\\system32"), otherwise the
program will not run.

For linux version, the GSL are included in the executables, so you can simply
put all files into one folder and just run the program by ./BEAM. 

[2. Input Format]-----------------------------------------

The user needs to create a file (by default, "data.txt") that contains the case-control 
genotype data to be used as an input to the program. Although our program is designed 
for biallelic SNPs, the user may use markers that take more than two alleles. However, 
BEAM assumes that all markers in the input file are of the same type (e.g., they are all 
k-allelic markers). The allele at a k-allelic marker should be indexed from 0 to (k-1). 
For example, without assuming Hardy-Wenberg Equilibrium (HWE), a genotype at a SNP locus 
takes 3 different values, which should be indexed as 0,1,2. It doesn't matter which allle
is coded by which integer, i.e., the user may code the wild type homozygote by 2, the 
mutant type homozygote by 1, and the heterozygote by 0. In addition, use negative numbers 
for missing alleles

The first line of the input file should contain the disease status of each individual.
You should use 1 to denote patients and 0 to denote controls. Alternatively, you may use 
a one-dimensional continuous trait value for each individual. BEAM will automatically 
assign those individuals with traits larger than the mean to cases, and the remaining 
individuals to controls. 

Starting from the second line of the input file, each line contains the genotype data
at one marker for all individuals, separated by ' ' (space). For diploid species, if 
you want to assume HWE in association mapping, then you should input two alleles at the
locus for each individual (the order of alleles doesn't matter). For example, for 2 
patients and 2 controls genotyped from 5 markers, if assuming HWE, the input file may 
look like:
1 1 1 1 0 0 0 0				<--disease status, two identical status per individual, 1: case, 0: control
1 0 0 0 1 1 0 0				<--marker 1, alleles are denoted by 0 and 1.
1 1 0 0 1 0 1 1				<--marker 2
0 0 0 1 0 0 0 1				<--marker 3
0 1 1 1 0 -1 -1				<--marker 4, "-1" denotes missing allele
1 1 0 1 1 0 1 0				<--marker 5

Since each column in the input file denote one individual (or a haploid),	the disease 
status for each individual (haploid) must match with the correponding column in the 
file. In the above example, when specifying two allels per individual, the user should 
input two identical disease status per individual	in the first line. 

The user may also provide the SNP ID and their chromosomal locations (in bps) at the 
begining of each line. The above example file then looks like:
ID 		Chr 	Pos 		1 1 1 1 0 0 0 0	
rs100102	chr1 	4924223 	1 0 0 0 1 1 0 0
rs291093	chr2	35981121	1 1 0 0 1 0 1 1
rs490232	chr9	6920101		0 0 0 1 0 0 0 1
rs093202	chrX	319101		0 1 1 1 0 -1 -1				
rs43229		chrY	103919		1 1 0 1 1 0 1 0				

IMPORTANT: 	if SNP ID and locations are provided in the data, make sure you include 
						in the first line "ID	Chr Pos " to label each column. This is 
						to ensure that the disease status of each individual matches to 
						the individual in the corresponding column. In addition, please use
						chr1,...,chr22,chrX, chrY to denote chromosomes. "Chr" is ok too.

Since BEAM is designed for detecting interaction between markers that are far apart,
by providing the marker location information, it helps BEAM to avoid being trapped
by local marker dependencies (due to LD rather than interaction effects). In addition, 
if the SNP ID and their locations are provided in the input file, you should turn on 
the correponding parameters in the file "parameter.txt". This can be done by setting 
both INC_SNP_ID and INC_SNP_POS to 1 (default).

Please see the "data.txt" provided online for an example of the input file.

[3. Program Parameters]-----------------------------------------

The running parameters of BEAM are specified in the "parameter.txt" file. These para-
meters include: input filename, outputfilename, priors, burnin-length, mcmc-length,
and thin, etc. The parameters "burnin", "mcmc", "thin" should be chosen according to 
the number of markers genotyped in the data. Denote the totoal number of markers by L, 
we suggest the following choice of parameters: 
burin = 10~100 L 
mcmc = L^2
thin = L

In addition, the user may let the program first search for some local modes of the joint 
posterior distribution, and then run MCMC from those modes. This strategy can be advantagous 
than starting from random points, espetically for detecting high-order interactions. To
do this, set INITIALTRYS to a positive number (e.g. 10 ~ 100), so that BEAM will first 
search for local modes in 10~100 trials. Set TRY_LENGTH to a number such as 10L for the
number of iterations in each trial (L is the number of markers). 

If during the MCMC run the algorithm finds a better local mode (with a larger likelihood) 
compared to the mode where BEAM starts with, you can let BEAM restart the chain from this 
better mode. To do this, set AUTORESTART to a number between 5~10, such that if the new 
local mode measured in log-likelihood is 5~10 larger than the initial mode, the chain will 
restart. To disable this function, set AUTORESTART to a large value, e.g., 1000000. Never 
set it to a negative number.

The user can let BEAM to automatically determine a set of parameters for a dataset, such as 
burnin-length, mcmc-length, thin, etc. By default, we set burnin = 100L if the MCMC chain 
starts from random configurations (i.e., when INITIALTRYS is 0). If the user let BEAM to 
search for some local modes first (by setting INITIALTRYS to a positive integer), we then 
set burnin = 10L. The length of each initial trial (TRY_LENGTH) is 20L by default. We further 
set mcmc = L * L, and thin = L. To use the automatic setting, the user must set the corresponding 
parameters in "parameter.txt" to 0 or a negative number. For example, let
BURNIN		0
MCMC		0
THIN		0
TRY_LENGTH	0

NOTE: Any positive integer will replace the default setting of BEAM.

BEAM will output markers/interactions with Bonferroni adjusted p-value smaller than a user-
specified value, which can be specified by "P_THRESHOLD" in "parameter.txt".

You can use BEAM to search for marginal associations only if let SINGLE_ONLY = 1.

You may specify the input file and the output file in "parameter.txt".

We encourage you to modify "parameter.txt" to run BEAM under different parameter settings. 
Please see "parameter.txt" for more details. 

[4. Command Line]-----------------------------------------

To run BEAM, type in the command line:

"BEAM [input output]"

BEAM will refer to "parameter.txt" for running parameters, so you can run the program
by simply typing "BEAM". The only option you can use in the command line is to 
specify the input file name and the output file name. Note that if you want to specify 
either of them, you must specify both of them.

[4. Output]-----------------------------------------

The main output file contains the estimated posterior probability of association for 
each marker, including both marginal and interactive association probabilities. The file 
also contains posterior estimates of the number of marginal associations and the size 
of interactions. In addition, we evaluate the p-value of each detected association using
the B-statistic introduced in the paper. Only significant p-values (<0.1 after the Bonferroni
correction) and associated markers are reported. More than one markers reported within 
a single parenthsis indicate interactions effects.

To check the performance of Markov chains, we output two more files: "lnp.txt" and
"posterior.txt". The former contains the trace of log-posterior probabilities (up to a 
normalizing constant) of the sampled parameters. The latter contains the summary 
of Markov chains. "posterior.txt" also contains the B-statistics and the estimated 
p-values for detected candidate markers and interactions.

[5. Credit]-----------------------------------------

This program is developed based on algorithms proposed in 

Yu Zhang, Jun Liu (2007) "Bayesian Inference of Epistatic Interactions in Case-Control
Studies", Nature Genetics, in press.

Please cite the paper if you used this program in your research for publication. 
The research was supported by an NIH R01 grant.
 
[6. Support]-----------------------------------------

All questions and comments should be directed to Yu Zhang at the Department of Statistics, 
The Pennsylvania State University, University Park, PA 16802.
Email: yuzhang@stat.psu.edu


#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work

$command = "smartpca.perl";
$command .= " -i example.geno ";
$command .= " -a example.snp ";
$command .= " -b example.ind " ;
$command .= " -k 2 ";
$command .= " -o example.pca ";
$command .= " -p example.plot ";
$command .= " -e example.eval ";
$command .= " -l example.log ";
$command .= " -m 5 ";
$command .= " -t 2 ";
$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "eigenstratQTL"; # or eigenstratQTL.big.perl for large data sets
$command .= " -i example.geno ";
$command .= " -j example.phenoQTL ";
$command .= " -p example.pca ";
$command .= " -l 1 ";
$command .= " -o example.chisqQTL ";
print("$command\n");
system("$command");

$command = "gc.perl example.chisqQTL example.chisqQTL.GC";
print("$command\n");
system("$command");

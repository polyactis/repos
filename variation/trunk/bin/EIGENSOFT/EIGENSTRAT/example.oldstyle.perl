#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}";

$command = "pca";
$command .= " -i example.geno ";
$command .= " -k 2 ";
$command .= " -o example.pca ";
$command .= " -e example.eval ";
$command .= " -l example.log ";
$command .= " -m 5 ";
$command .= " -t 2 ";
$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "eigenstrat";
$command .= " -i example.geno ";
$command .= " -j example.pheno ";
$command .= " -p example.pca ";
$command .= " -l 1 ";
$command .= " -o example.chisq ";
print("$command\n");
system("$command");

$command = "gc.perl example.chisq example.chisq.GC";
print("$command\n");
system("$command");

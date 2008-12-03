#!/usr/bin/perl

$command = "../bin/smartpca";
$command .= " -p par.example >example.log";
print("$command\n");
system("$command");

$command = "../bin/ploteig";
$command .= " -i example.evec ";
$command .= " -c 1:2 ";
$command .= " -p Case:Control ";
$command .= " -x ";
$command .= " -o example.plot.xtxt "; # must end in .xtxt
print("$command\n");
system("$command");

$command = "../bin/evec2pca.perl 2 example.evec example.ind example.pca";
print("$command\n");
system("$command");

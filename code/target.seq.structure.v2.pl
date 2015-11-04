#!/usr/bin/perl -w
# Name:				target.seq.structure.v2.pl (study the sequence structure for target regions)
# Organization:		Case Western Reserve University
# Date Updated:  	9/8/2010 (Shuying made minor changes on comments) 
#                   9/14/2010 Shuying modified Stephen's original CpGI.structure.pl to make it 
#                   possible to check the structure of all general target region. 
# Change Log:		1) Using regular expression for finding cg's
# Assumptions:		This perl script assumes that you name your reference files as "chrN.fa". 
#					For example, one could set reference as "chr3.fa". Furthermore, it assumes that
#					for each CpG island identified within the query file (second), each chromosome 
#					parameter) found will have a corresponding reference FASTA file in the 
#					directory path specified.
#Purpose:			This Perl script is designed to analyze target region sequence information and hg 
#					FASTA references for the information is is described in "output"

# Input:			perl target.seq.structure.v2.pl </path/to/referenceFASTA/directory> <path/to/CpG/location/file> <headers T|F> <out name>
# Example:			perl target.seq.structure.v2.pl /home/projects/data/reference/hg18/ /home/projects/data/non.seq.data/demoX.bait.wig  F
#						out.file.name
# Output (20 columns) :
#   1) chr  start  end  length
#   2) A, C, G, T, C+G, CG, nonCG-C,
#   3) percentage of A, C, G, T, C+G, CG, nonCG-C,
#   4) lower_count		%lower_count
use warnings;
use strict;


#file usage, checks for valid number of input parameters
my $usage = qq(Usage:  perl structure.pl </path/to/reference/directory> <path/to/CpG/location/file> <headers T|F>);

die($usage) if (@ARGV < 4);

#create the file handlers and open
open(ISLANDS, $ARGV[1]) || die("the path is incorrect") or die $!; #$ARGV[1]
open(OUTPUT, ">".$ARGV[3]) or die $!;	 #the date will be attached to the prefix of the output

#output header (column names)
print OUTPUT "chr		start	end		len		a	c	g	t	c+g cg 	ncgc pa		pc		pg		pt		pc+g	pcgc	pncgc	lower_count		%lower_count\n";

#global headers
my $island; my @islandSplit;
my $numIslands = 0; my $chr; my $start; my $end;
my $length; my $fHandle; my $bool; my $linesToSkip;
my $temp; my $seq; my $cursor; my $char; my $sub;
my $a; my $c; my $g; my $t; my $cg; my $cpg; my $noncgc;
my $pa; my $pc; my $pg; my $pt; my $pcg; my $pcpc; my $pnoncgc; my $pn; my $obsexp; my $islands; my $cgc; my $pcgc;
my $lower_count; my $plower_count; sub newChr;
my $lastChr = ""; my $reference; my $pcpg;

#skip the column labels line
if($ARGV[2] eq "T"){
	$islands = <ISLANDS>;
}

while($island = <ISLANDS>){ #while = for each CpG island, do the following:
	$numIslands++;
	$a = 0; $c = 0; $g = 0; $t = 0; $cg = 0; $cgc = 0; $noncgc = 0; $lower_count = 0;
	chomp $island;
	@islandSplit = split("\t", $island);
	$chr = $islandSplit[0];
	$start = $islandSplit[1];
	$end = $islandSplit[2];	
	$length = $end - $start + 1;

	#if a new chr is found, generate the new reference sequence
	if($chr ne $lastChr){
		$reference = newChr($chr);
		print("Found new chromosome $chr\n");
	}	
	$lastChr = $chr;
	#generate the new target sequence
	$seq = substr($reference, $start - 1, $length);
	if($start == 3320900){
		print("seq=> $seq\n");
	}
	$lower_count = $lower_count + ($seq =~ tr/a/A/);
	$lower_count = $lower_count + ($seq =~ tr/c/C/);
	$lower_count = $lower_count + ($seq =~ tr/g/G/);
	$lower_count = $lower_count + ($seq =~ tr/t/T/);
	$a = ($seq =~ tr/A/A/);  #count the bases
	$c = ($seq =~ tr/C/C/); 
	$g = ($seq =~ tr/G/G/); 
	$t = ($seq =~ tr/T/T/);
	$cg = $seq =~ s/CG/CG/gi;
	
	#calculate the percentages
	if($length != 0 and $c != 0){
		$pa = $a / $length * 100; $pc = $c / $length * 100; $pg = $g / $length * 100; $pt = $t / $length * 100; $noncgc = $c - $cg;
		$cpg = $c + $g; $pcpg = $cpg / $length * 100; $pcg = $cg / $c * 100; $pnoncgc = $noncgc / $c * 100;
		$obsexp = $cg * $length / ($c * $g); $plower_count = $lower_count / $length * 100;

	}
	
	printf OUTPUT "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\n", $chr, $start, $end, $length, $a, $c, $g, $t, $cpg,$cg, $noncgc, $pa, $pc, $pg, $pt, $pcpg, $pcg, $pnoncgc,  $lower_count, $plower_count;	
}

close ISLANDS; close OUTPUT;

sub newChr(){
	my $fHandle;
	my $sequence;
	my $line;
	if($ARGV[0] =~ m/\/$/){
		$fHandle = $ARGV[0].$_[0].".fa";
	}
	else{
		$fHandle = $ARGV[0]."/".$_[0].".fa";
	}
	open(REBUILD, $fHandle) or die("Failed to open the reference FASTA file\n");
	$line = <REBUILD>;
	while($line = <REBUILD>){
		chomp($line);
		$sequence = $sequence.$line;
	}
	close REBUILD;
	return $sequence;
}

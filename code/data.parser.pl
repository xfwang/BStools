# Remove this "#!/usr/bin/perl -w" in the fist line  to avoid some unnecessary error message "Args must match #1". 
# data.parser.pl 
use strict; use warnings;

# perl data.parser.pl *forw.txt output_dir prefix


my $date1 = `date`;
print "Start separating the base file: $date1\n";
open (bOUT, '>', "$ARGV[1]/$ARGV[2].base") or die "unable to write .base $!\n";

open(POS, "<$ARGV[0]") or die "error reading file";

while (my $base_line =<POS>)
 {
 	chomp($base_line); 
 	my @line_split=split("\t", $base_line); 
 	   my $C_col = $line_split[3]; 
       my @C_array = split(":", $C_col); # to separate the CG-type from coverage by colon
       my $C_type = $C_array[0]; 
   	   my $total_CT = $C_array[1];
       my $start = $line_split[1] + 1; # offset the start position by one (1-based offset)
 	
        	print bOUT "$line_split[0]\t$start\t$C_type\t$total_CT\t$line_split[4]\n";
 }
close bOUT;
my $date2 = `date`;
print "Finished separating the base file: $date2\n";



my $date3 = `date`;
print "Start separating nonCG and CG sites: $date3\n";

open(BASE, "<$ARGV[1]/$ARGV[2].base") or die "error reading file"; # base input file
open(BNCGC, ">$ARGV[1]/$ARGV[2].nonCG.base") or die "error writing to file"; # BNCGC -> base non CGc
open(BNCGCwC, ">$ARGV[1]/$ARGV[2].nonCG.covered.base") or die "error writing to file"; # BNCGCwC -> base non CGc w/ coverage
open(BCGC, ">$ARGV[1]/$ARGV[2].CG.base") or die "error writing to file"; # BCGC -> base CGc

while (my $base_line = <BASE>) {

	chomp($base_line);
	my @base_split = split("\t", $base_line);
	
	if ($base_split[2] ne "CG") {
		print BNCGC "$base_split[0]\n";
		
		if ($base_split[3] > 0) {
  			print BNCGCwC "$base_split[0]\t$base_split[4]\n";
  		}
  	}

	else {
        print BCGC "$base_split[0]\t$base_split[1]\t$base_split[3]\t$base_split[4]\n";

	}
}#end while

my $date4 = `date`;
print "Finished separating nonCG and CG sites: $date4\n";


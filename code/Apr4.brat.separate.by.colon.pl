
#! This perl file is used to separate the brat acgt -B output file. 
# and also convert the position from 0-based offset to 1-based offset
# The input file is a file with the following lines
#  chrX	0	0	CHH:0	0	+
#  chrX	4	4	CHH:0	0	+

print "--------------------------------------------- \n"; 
print "The 7 columns in the output files are: \n"; 
print "chr,  start (1-based), end (1-based), C_type, total_reads, methy_level, and strand \n"; 
print "--------------------------------------------- \n"; 
print "The usage of this perl code file is:  perl  Apr4.brat.separate.by.colon.pl  B.acgt_forw.txt  output \n"; 

use strict; use warnings;

open (IN, $ARGV[0])  or die("Failed to open the input file at $ARGV[0] \n"); 
open (OUT, ">$ARGV[1]") or die("Can not open output file \n"); 


my $j = 1; 
while (my $line1 = <IN>) 
{ 
   if ( $j%50000000 ==0 ) 
   { my $now = localtime time;
      print(" when j is: $j, time is:, $now. \n") ; 
   }
   $j++; 

   chomp($line1); 
   my @line_split=split("\t", $line1); 
   my $C_col = $line_split[3]; 
   my @C_array = split(":", $C_col);
   my $C_type = $C_array[0]; 
   my $total_CT = $C_array[1];
   my $start = $line_split[1] + 1; 
   my $end = $line_split[2] + 1 ; 

   print OUT "$line_split[0]\t$start\t$end\t$C_type\t$total_CT\t$line_split[4]\t$line_split[5]\n";  

} 

close IN; 
close OUT; 

# Remove this "#!/usr/bin/perl -w" in the fist line  to avoid some unnecessary error message "Args must match #1".  
# bisulfate.extract.reads.v2.pl
use strict; use warnings;

############################################################################################################
#	Purpose: This code's purpose is to separate input data into good and bad BS rate files based on the condition
#			that a good read has median BSrate > 0.99 and a bad read has median BSrate <= 0.6. The uses can change 
#           these settings in the GLOBAL PARAMETER secion of the MethyQA.pl and partial.MethyQA.pl 
#	
#	File: 	bisulfate.extract.reads.pl
#	Usage: 	perl bisulfate.test.pl INPUT High_BS_Rate Low_BS_Rate High_Cov_Rate Low_Cov_Rate
#
#	Input:	Input file needs to be a summary table generated for target dataset
#			Example located at:
#			cd /home/ajn21/bisulfate/output/target/for/forward.target.summary.txt
#	
#	Note0: 	If you view the above file you'll notice that the first 3 columns are technically
#			one column, as they serve as the ID for the summary data. For example:
#			          ID               Min 1stQu Median Mean 3rdQu Max
#			chr1 112287034 112287223	0.8	1	1	0.9931	1	1
#			chr1 119392928 119393132	1	1	1	1	1	1
#			chr1 145538469 145538531	1	1	1	1	1	1
#			chr1 145847645 145847830	1	1	1	1	1	1
#
#			The column names aren't listed in the *summary.txt files but the 6 numbers after ID result from the 
#			summary() function in R. Therefore, col0 is the ID, col1 is the Min, col2 is the 1stQu, col3 is the Median,
#			col4 is the Mean, col5 is the 3rdQu, and col6 is the Max. In this case the code separates data based on the median
#			value in column 3
#
#	
#	Output: The output will be 2 files with same format as input except with separated median values
#
#	Example:perl /home/ajn21/bisulfate/bisulfate.extract.reads.pl /home/ajn21/bisulfate/output/target/rev/RG1/reverse.target.summary.txt 
#				 0.9  0.1  0.8  0.2
# 
############################################################################################################

open(FILE, "<$ARGV[0]") or die "error reading file";
open(highBS, ">$ARGV[1]") or die "error opening $ARGV[1]";
open(lowBS, ">$ARGV[2]") or die "error opening $ARGV[2]";
open(highP, ">$ARGV[3]") or die "error opening $ARGV[3]";
open(lowP, ">$ARGV[4]") or die "error opening $ARGV[4]";


while (my $line = <FILE>) {

	chomp($line);
	my @split = split("\t", $line);
	
	if ($split[8] >= $ARGV[5]) { # median >= cutoff
	print highBS "$line\n";
	}
	if ($split[8] <  $ARGV[6]) { # median <= cutoff
	print lowBS "$line\n";
	}
	if ($split[5] >= $ARGV[7]) { # coverage > cutoff%
	print highP "$line\n";
	}
	if ($split[5] <  $ARGV[8]) { # coverage <  cutoff%
	print lowP "$line\n";
	}
		
}#end while
close FILE;
close highBS;
close lowBS;
close highP;
close lowP;

exit;
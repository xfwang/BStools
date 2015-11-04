# Remove this "#!/usr/bin/perl -w" in the fist line  to avoid some unnecessary error message "Args must match #1". 
# BQA.data.parser.v4.pl 
use strict; use warnings;

# Note, originally, the first line of this code is "#!/usr/bin/perl -w". Since sometimes
# we got error message " Args must match #1". I remove it. 

#############################################################################
#	Purpose: This code's purpose is to separate input data with
#			multiple chr values into separate files with only
#			one chr value in each output file. After separation,
#			a call to a subroutine is made to match each base position
#			to each target region for each respective chromosome. Then
#			the output files are parsed to separate nonCGc sites.
#	
#	v3:		Version 3 of the BQA data parser combines the Apr4.brat.separate.by.colon.pl
#			with version 2 of the data parser. 
#   v4:     To speed up this parsing program when dealing with large 19GB ACGT file, I (ss)
#           add if ( $line_split[0] eq $ARGV[2]) {  } to only deal with one specific chromsome
#           This change make the program 15 minutes faster (51 vs 36 min) 
#	
#	File: 	BQA.data.parser.v3.pl
#	Usage: 	BQA.data.parser.v3.pl Target_INPUT Base_INPUT chr Output_Directory prefix 
#	
#	Input:	
#			Target_INPUT - a file with target regions listed
#			  chr start end 
#			chr1	496	588	
#			chr1	648	719	
#			chr1	720	795	
#			chr1	796	919	
#
#			Set the target input as "F" if the user does not have a target region file
#
#			Base_INPUT - the acgt-count output file
#			chr1	50006	50006	CHH:0	0	+
#			chr1	50007	50007	CHH:0	0	+
#			chr1	50014	50014	CHH:0	0	+
#			chr1	50015	50015	CHH:0	0	+
#
#
#	Output: There will be multiple outputs from this script. The *.base files will
#			represent individual chromosomes separated from the acgt-count file. The *.out
#			files will represent the base positions mapped to target regions and will be a 
#			combination of the *.base and *.target files. The *.nonCG files will represent
#			all nonCGc sites separated from CGc sites, and *.nonCGc.coverage files represent
#			only nonCGc sites with coverage greater than zero.
#
#	Note 1:	Part2 and Part3 (see below) can be done using mawk:			
#			date; mawk '{outFile=$1; print $0 >outFile}' ../B.Mar31.BS1.acgt_forw.1base.Offset ; date
#
#############################################################################

if ($ARGV[0] ne 'F') {
########################### sorting the target input file ########################

##########
# PART1:
##########
# The following step is used to sort the target input file.
my $date = `date`;
print "Sorting the target file $ARGV[0]: $date\n"; 
`sort -t \$'\t'  -k 1,1  -k 2n,2  $ARGV[0] > $ARGV[0].sorted`; # first sort by chromosome, then sort by start position.
my $date1 = `date`;
print "Finished sorting the target file: $date1\n";

#############################################################################

} # end of target flag conditional


########################### separation of chr values ########################

if ($ARGV[0] ne 'F') {

##########
# PART2:
##########
# The following step is used to get all data for each chr from the target file.

# Using the Target_Input, Nov18.CCGG.all.40to220bp.w.status.minus1.TAB (see notes above)
# Part1 takes approximately 4 minutes to complete

my $date2 = `date`;
print "Start separating the target file: $date2\n";
my $curr_chrom_reg = "$ARGV[2]"; # will serve as my marker
open (rOUT, '>', "$ARGV[3]/$ARGV[4].$curr_chrom_reg.target") or die "unable to write $curr_chrom_reg.target $!\n";

open(REG, "<$ARGV[0].sorted") or die "error reading file";
while (<REG>)
{
   my ($chrom_reg) = split;  # $chrom_reg is the first column
   if ($chrom_reg eq $curr_chrom_reg) 	
   {								
    print rOUT;
   }
}
close rOUT;
my $date3 = `date`;
print "Finished separating the target file: $date3\n";

} # end of target flag conditional

##########
# PART3:
##########
# The following step is used to get all data for each chr from the base file.

# Using the Base_Input, B.Mar31.BS1.acgt_forw.1base.Offset, (see notes above)
# Part2 takes approx. 11 minutes to complete
my $date4 = `date`;
print "Start separating the base file: $date4\n";
my $curr_chrom_base = "$ARGV[2]";
open (bOUT, '>', "$ARGV[3]/$ARGV[4].$curr_chrom_base.base") or die "unable to write $curr_chrom_base.base $!\n";

open(POS, "<$ARGV[1]") or die "error reading file";

while (my $base_line =<POS>)
 {
 	chomp($base_line); 
 	my @line_split=split("\t", $base_line); 
 	if ( $line_split[0] eq $ARGV[2]) 
 	{  
 	   my $C_col = $line_split[3]; 
       my @C_array = split(":", $C_col); # to separate the CG-type from coverage by colon
       my $C_type = $C_array[0]; 
   	   my $total_CT = $C_array[1];
       my $start = $line_split[1] + 1; # offset the start position by one (1-based offset)
       my $end = $line_split[2] + 1 ;  # offset the end position by one (1-based offset)
 	
    	my ($chrom_base) = $line_split[0];  # $chrom_base is the first column
    	if ($chrom_base eq $curr_chrom_base) # if a new chromosome is discovered, a new file is written
    	{
        	print bOUT "$line_split[0]\t$start\t$end\t$C_type\t$total_CT\t$line_split[4]\t$line_split[5]\n";
    	}
    }
 }
close bOUT;
my $date5 = `date`;
print "Finished separating the base file: $date5\n";


########################  identify then call to subroutine  ####################

if ($ARGV[0] ne 'F') {

##########
# PART4:
##########
# This step is used to identify if file was created in parts 1 and 2
# and if it was, then run target and base data for each chr to subroutine
# to find base position in target region.


my $date6 = `date`;
print "Start finding base in target region: $date6\n";

my $key = "$ARGV[2]";
my $target_file = "$ARGV[4].$key.target";
my $base_file = "$ARGV[4].$key.base";
	if (-e $target_file) { #identifies if the file was created in the above separation process
		if (-e $base_file) {
			my $arg0 = "$ARGV[4].$key.target"; #temp variables used for call to subroutine
			my $arg1 = "$ARGV[4].$key.base";
			my $arg2 = "$ARGV[4].$key.out";
			posTarget($arg0, $arg1, $arg2); #call to subroutine
			my @timeData = localtime(time);
			print "$timeData[2]:$timeData[1]:$timeData[0]\n"; # so I can view the time after each cycle
			print "finished processing $key\n";				  # Hour:Minute:Second
			print "-------------------------\n";
		}
	}
	

} # end of target flag conditional

########################  Filter nonCGc and nonCGc w/ coverage  ####################

##########
# PART5:
##########
open(BASE, "<$ARGV[3]/$ARGV[4].$ARGV[2].base") or die "error reading file"; # base input file
open(BNCGC, ">$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.base") or die "error writing to file"; # BNCGC -> base non CGc
open(BNCGCwC, ">$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.covered.base") or die "error writing to file"; # BNCGCwC -> base non CGc w/ coverage
if ($ARGV[0] ne 'F') {
open(TARG, "<$ARGV[3]/$ARGV[4].$ARGV[2].out") or die "error reading file"; # target (paired w/ base pos) input file
open(TNCGC, ">$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.out") or die "error writing to file"; # TNCGC -> target non CGc
open(TNCGCwC, ">$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.covered.out") or die "error writing to file"; # TNCGCwC -> target non CGc w/ coverage
} # end of target flag conditional

while (my $base_line = <BASE>) {

	chomp($base_line);
	my @base_split = split("\t", $base_line);
	
	if ($base_split[3] ne "CG") {
		print BNCGC "$base_split[5]\n";
		
		if ($base_split[4] > 0) {
  			print BNCGCwC "$base_split[0]\t$base_split[5]\n";
  		}
  	}
}#end while


if ($ARGV[0] ne 'F') {
while (my $targ_line = <TARG>) {

	chomp($targ_line);
	my @ts = split("\t", $targ_line);

	if ($ts[5] ne "CG") {
	#if ($ts[5] eq "CG") {
		print TNCGC "$ts[0]\t$ts[1]\t$ts[2]\t$ts[6]\t$ts[7]\n";
		
		if ($ts[6] > 0) {
  			print TNCGCwC "$ts[0]\t$ts[1]\t$ts[2]\t$ts[6]\t$ts[7]\n"; 								
  		}
  	}  		
}#end while

} # end of target flag conditional

close BASE;
close BNCGC;
close BNCGCwC;

if ($ARGV[0] ne 'F') {
close TARG;
close TNCGC;
close TNCGCwC;



##########
# PART6:
##########

open(COV, "<$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.covered.out") or die "error reading file";
my %Chrm;
my (@Cov,@line);
while (my $cov = <COV>) {
	chomp $cov;
	@Cov=split("\t",$cov);
	$Chrm{$Cov[1]} = undef; 
}


open(FULL, "<$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.out") or die "error reading file";
open(MATCH, ">$ARGV[3]/$ARGV[4].$ARGV[2].nonCG.match") or die "error reading file";

while (my $line=<FULL>) {
	
	chomp $line;
	@line=split("\t",$line);
	
	if (exists $Chrm{$line[1]} ) {
		print MATCH "$line\n";
	}	
}

close COV;
close FULL;
close MATCH;

} #end of target flag conditional

##############################  end of main process  ############################

########################  start of subroutine  ##################################

sub posTarget {

my($arg0, $arg1, $arg2) = @_;
open(TARGET, "<$ARGV[3]/$arg0") or die "error reading $arg0";
open(BASE, "<$ARGV[3]/$arg1") or die "error reading $arg1";
open(OUT, ">$ARGV[3]/$arg2") or die "error reading $arg2";

  my @base_arr;
   my @target_arr;
   my $base_line; 
   my $target_line;

   # while(my $target_line=<TARGET>)
   while($target_line=<TARGET>)
   {	
        chomp($target_line);
	while( $base_line=<BASE>)
		{
			if(!$target_line)
			{
			    last; #This is done because once we reach the end of the target file the loop should be terminated
			}
                        chomp($base_line);
			@target_arr=split("\t",$target_line);
			@base_arr=split("\t",$base_line);
			
			if($base_arr[1] > $target_arr[2])
			{
			    $target_line=<TARGET>; 
                            if(!$target_line)
		       	    {
				last;
			    }
                            # Without the above three lines, I will get the following error message while running the next "chomp($target_line)"
                            # Use of uninitialized value in scalar chomp at /home/projects/perl.code/ss.perl/June25.get.base.from.target.diff.chr.pl line 114, <TARGET> line 2.
                            chomp($target_line);  

			    redo;
			}
			elsif($base_arr[1] >= $target_arr[1] && $base_arr[1] <= $target_arr[2])
			{   # shift(@base_arr); # this save the one with the first column (chr) removed. 
                            # print OUT join ("\t", $target_arr[0], $target_arr[1], $target_arr[2], @base_arr), "\n"; 
                            # The above line will give extra blank line, so I use the following one:
                            # print OUT join ("\t", $target_arr[0], $target_arr[1], $target_arr[2], $base_arr[1], $base_arr[2], $base_arr[3], $base_arr[4], $base_arr[5], $base_arr[6]), "\n"; 
                            print OUT "$target_arr[0]\t$target_arr[1]\t$target_arr[2]\t$base_arr[1]\t$base_arr[2]\t$base_arr[3]\t$base_arr[4]\t$base_arr[5]\t$base_arr[6]\n"; 
			}                                           
		}
    }
close BASE;
close TARGET;
close OUT;
}# end of subroutine

###########################end of subroutine##################################
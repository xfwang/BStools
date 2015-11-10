# Remove this "#!/usr/bin/perl -w" in the fist line  to avoid some unnecessary error message "Args must match #1".  
# MethyQA.pl
use strict; use warnings;

###########################################################################################
#	Purpose: This script is the MethyQA pipeline, conducting steps 1-5 quality assessment as mentioned in the MethyQA manuscript
#
#
#	File: 	 BStools.pl
#	Usage: 	 perl BStools.pl -i <FASTQ_input> -t <TARGET_input> -d <MethyQA_directory>  -p <prefix> -R <reference_directory> -r <reference_name> -Y <python_pathway>
#
#
###########################################################################################
 
###################
# Global Variables
###################

my $old = select(STDOUT);
$| = 1;
select($old); 
my $numArgs = $#ARGV + 1;
use Cwd;
use Getopt::Std;
my %Options;
my $ok = getopts('i:t:d:p:R:r:f:a:A:T:N:n:q:m:H:F:B:b:L:l:u:v:Y:Q:M:K:J:X:C', \%Options);

# 14 args is the minimun number of options for the pipeline to work correctly

if($numArgs < 12 || !$ok || !exists($Options{i}) || !exists($Options{t}) || !exists($Options{d})  
||  !exists($Options{p}) || !exists($Options{R}) || !exists($Options{r})   ) 
{ 
  printUsage();
  exit;
}



############# Setting Global Parameters #####################

my $input = $Options{i};  # fastq input file
unless (-e $input) {
 print "Fastq File Doesn't Exist!\n\n";  exit 
 } 


my $target = $Options{t}; # target region input file
if ($target ne "F") 
  {
      unless (-e $target) 
		  {  print "Target File Doesn't Exist!\n\n";   exit } 
  }

my $directory = $Options{d}; # the path to MethyQA directory
my $output = cwd(); # output directory


# my $chr = $Options{c};    # chromosome number (e.g. chr2, chr5, chrX, chrY)
my $prefix = $Options{p}; # prefix written to output filenames
 

my $code_directory = $directory . '/code/'; # the directory where the code files for the pipeline are located

print "the dirctory is $directory \n";
print "the code dirctory is $code_directory \n"; 
print "the output dirctory is $output \n"; 

my $reference_directory = $Options{R}; # the reference directory, used for step 5
my $reference = $Options{r}; # user specifies the reference genome, used in step 3


my $python_version="";

my $fastqc_dir = ""; # path to fastqc
if ( !exists($Options{Q}) )
  { $fastqc_dir = $directory . '/resources/FastQC/./fastqc' ; } # use the fastqc compiled in MethyQA 
else 
  { 
	  	unless (-e $Options{Q}) 
		  {  print "Path to Fastqc Doesn't Exist!\n\n";   exit } 
	    $fastqc_dir = $Options{Q} ; 
}
 



my $brat_trim_dir = ""; # path to brat_trim
if ( !exists($Options{M}) )
  { $brat_trim_dir =  $directory . '/resources/brat/./trim'; } # use the brat-trim compiled in MethyQA 
else 
  {  
	     unless (-e $Options{M}) 
		  {  print "Path to BRAT_trim Doesn't Exist!\n\n";   exit } 
          $brat_trim_dir = $Options{M}
  }

my $brat_large_dir = ""; # path tobrat_large
if ( !exists($Options{K}) )
  { $brat_large_dir =  $directory . '/resources/brat/./brat-large'; } # use the brat-large compiled in MethyQA 
else 
  {  
	     unless (-e $Options{K}) 
		  {  print "Path to BRAT_large Doesn't Exist!\n\n";   exit } 
          $brat_large_dir = $Options{K}
  }


my $brat_acgt_dir  = ""; # path tobrat_large
if ( !exists($Options{J}) )
  { $brat_acgt_dir  =  $directory . '/resources/brat/./acgt-count'; } # use the brat-large compiled in MethyQA 
else 
  {  
	     unless (-e $Options{J}) 
		  {  print "Path to BRAT_acgt_count Doesn't Exist!\n\n";   exit } 
          $brat_acgt_dir = $Options{J}
  }


my $fastx_dir = ""; # path to fastx_clipper
if ( !exists($Options{X}) )
  { $fastx_dir =  $directory . '/resources/fastx/bin/./fastx_clipper';
    my $fastx_quality_trim_dir =  $directory . '/resources/fastx/bin/./fastq_quality_trimmer';
    my $fastx_trim_dir =  $directory . '/resources/fastx/bin/./fastX_trimmer';

  }  # use the fastx compiled in MethyQA 
else 
  {  
	     unless (-e $Options{X}) 
		  {  print "Path to Fastx_clipper Doesn't Exist!\n\n";   exit } 
          $fastx_dir = $Options{X}
  }




my $fastq_format = "";  
# If it is 'sanger' format fastq data, when trimming reads we set "-L 33";                               
# if it is 'illumina', we set "-L 64". The default is sanger. 
if ( !exists($Options{f}) )
  { $fastq_format = "sanger" }
elsif ( $Options{f} eq "sanger")
  { $fastq_format = "sanger" }
elsif ($Options{f} eq "illumina")
  { $fastq_format = "illumina" }
else
  { print "\n\nERROR: -f option should be set to 'sanger' or 'illumina', default is 'sanger'\n\n"; exit }
  

my $adapter_trim = "";
# If user specifies fastx or cutadapt as their desired adapter trimmer, then the user must provide
# the path of the programs using option -X for fastx and -x for cutadapt
if ( !exists($Options{a}) || $Options{a} eq "no" )
  { $adapter_trim = "no" }
elsif ( $Options{a} eq "fastx")
  { $adapter_trim = "fastx" }
elsif ($Options{a} eq "cutadapt")
  { 
	  $adapter_trim = "cutadapt";
      if (!exists($Options{Y})) 
		 {
		     print "\n\nERRPR: if -a is set as cutadapt, -Y <path_python> has to be specified in command line\n\n";
             printUsage();
             exit;
         }
	  else {$python_version= $Options{Y};} # user specifies how they call Python in their system
    }
else
  { print "\n\nERROR: -a option should be set to 'no' (default), 'fastx', or 'cutadapt' for no adapter trimming\n\n"; exit }


my $cutadapt_dir = ""; # path to fastx_clipper
if ( !exists($Options{C}) )
  { $cutadapt_dir =  $python_version . ' '. $directory . '/resources/cutadapt/bin/./cutadapt';} # use the fastx compiled in MethyQA 
else 
  {  
	     unless (-e $Options{C}) 
		  {  print "Path to Cutadapt Doesn't Exist!\n\n";   exit } 
          $cutadapt_dir = $python_version . ' '. $Options{C} 
  }


my $adapter_sequence = "";
### fastq format should be either "sanger" or "illumina"; default is "sanger"
if ( !exists($Options{A}) )
  { $adapter_sequence = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG" }
else
  { $adapter_sequence = $Options{A} }
  


## Quality Trimming:

my $trim_Flag = ""; 
# set as 'brat' if using BRAT trim function
# set as 'fix' if using fixed trimming procedure
# set as 'no' if the user does not want to do any quality trimming
if ( !exists($Options{T}) || $Options{T} eq "brat" )
  { $trim_Flag = "brat" }
elsif ( $Options{T} eq "fix")
  { $trim_Flag = "fix" }
elsif ( $Options{T} eq "retrieve")
  { $trim_Flag = "retrieve"}
elsif ( $Options{T} eq "no")
  { $trim_Flag = "no" }
else
  { print "\n\nERROR: -T option should be set to 'fix', 'retrieve', 'no', or 'brat' (default) for the type of trimming\n\n"; exit }
  
  
my $fixed_trim_5;	
# for fixed trimming, specifies the number of bases to trim at 5' end
if ( !exists($Options{N}) )
  { $fixed_trim_5 = 5 }
elsif ( $Options{N} >= 0 )
  { $fixed_trim_5 = $Options{N} }
else
  { print "\n\nERROR: -N option should be set to integer >=0, default is 5\n\n"; exit }
  
  
my $fixed_trim_3;	
# for fixed trimming, specifies number of bases to trim at 3' end
if ( !exists($Options{n}) )
  { $fixed_trim_3 = 10 }
elsif ( $Options{n} >= 0 )
  { $fixed_trim_3 = $Options{n} }
else
  { print "\n\nERROR: -n option should be set to integer >=0, default is 10\n\n"; exit }
  
my $brat_trim_qc;	
# for brat trimming, specifies the threshold for quality score
if ( !exists($Options{q}) )
  { $brat_trim_qc = 10 }
elsif ( $Options{q} >= 0 )
  { $brat_trim_qc = $Options{q} }
else
  { print "\n\nERROR: -q option should be set to integer >=0, default is 10\n\n"; exit }

my $brat_trim_N;	
# for brat trimming, specifies the allowed number of Ns
if ( !exists($Options{m}) )
  { $brat_trim_N = 2 }
elsif ( $Options{m} >= 0 )
  { $brat_trim_N = $Options{m} }
else
  { print "\n\nERROR: -m option should be set to integer >=0, default is 2\n\n"; exit }

my $brat_mismatch;  
# for brat-large, specifies the threshold for mismatch
if ( !exists($Options{H}) )
  { $brat_mismatch = 2 }
elsif ( $Options{H} >= 0 )
  { $brat_mismatch = $Options{H} }
else
  { print "\n\nERROR: -H option should be set to integer >=0, default is 2\n\n"; exit }


my $seed;  
# for brat-large, specifies the length of first F bases allow NO mismatches
if ( !exists($Options{F}) )
  { $seed = 24 }
elsif ( $Options{F} >= 0 )
  { $seed = $Options{F} }
else
  { print "\n\nERROR: -F option should be set to integer >=0, default is 2\n\n"; exit }


### Step 5 Variables:

my $high_bisulfite; 
# cutoff for high bisulfite rate for use in step 5
if ( !exists($Options{B}) )
  { $high_bisulfite = 0.99 }
elsif ( $Options{B} >= 0.0 && $Options{B} <= 1.0 )
  { $high_bisulfite = $Options{B} }
else
  { print "\n\nERROR: -B option should be set to float (0.0 <= B <= 1.0), default is 0.99\n\n"; exit }
  
  
my $low_bisulfite;  
# cutoff for low bisulfite rate for use in step 5
if ( !exists($Options{b}) )
  { $low_bisulfite  = 0.6 }
elsif ( $Options{b} >= 0.0 && $Options{b} <= 1.0 )
  { $low_bisulfite = $Options{b} }
else
  { print "\n\nERROR: -b option should be set to float (0.0 <= B <= 1.0), default is 0.99\n\n"; exit }
  
  
my $high_coverage;  
# cutoff for high coverage for use in step 5. Here, 0.5 means that at least
# 50% of CG sites with at least 1X coverage in a target region 
if ( !exists($Options{L}) )
  { $high_coverage  = 0.5 }
elsif ( $Options{L} >= 0.0 && $Options{L} <= 1.0 )
  { $high_coverage = $Options{L} }
else
  { print "\n\nERROR: -L option should be set to float (0.0 <= B <= 1.0), default is 0.5\n\n"; exit }

  
my $low_coverage; 
# cutoff for low coverage for use in step 5. Here, 0.1 means that at most
# 10% of CG sites with at least 1X coverage in a target region 
if ( !exists($Options{l}) )
  { $low_coverage   = 0.1 }
elsif ( $Options{l} >= 0.0 && $Options{l} <= 1.0 )
  { $low_coverage = $Options{l} }
else
  { print "\n\nERROR: -l option should be set to float (0.0 <= B <= 1.0), default is 0.1\n\n"; exit }

  
my $bisulfite_Flag = "";
# flag to indicate whether to make a boxplot to compare the sequence structure of 
# target regions with high or low bisulte rates.
if ( !exists($Options{u}) || $Options{u} eq "TRUE" )
  { $bisulfite_Flag = "TRUE" }
elsif ($Options{u} eq "FALSE") 
  { $bisulfite_Flag = "FALSE" }
else 
  { print "\n\nERROR: -u option should be logic value, either 'TRUE'(default) or 'FALSE' \n\n"; exit }

    
my $coverage_Flag  = "";  
# flag to indicate whether to make a boxplot to compare the sequence structure of 
# target regions with high or low coverage.
if ( !exists($Options{v}) || $Options{v} eq "TRUE" )
  { $coverage_Flag = "TRUE" }
elsif ($Options{v} eq "FALSE") 
  { $coverage_Flag = "FALSE" }
else 
  { print "\n\nERROR: -v option should be logic value, either 'TRUE'(default) or 'FALSE' \n\n"; exit }



############
# Step 1:
############
# Time: 
# Purpose: run fastqc to do quality assessment
system("date;");
print "------------------------------------------------------------------------\n";
print "Starting Step 1: FastQC\n";
print "------------------------------------------------------------------------\n\n";
system("$fastqc_dir $input -o $output");



############
# Step 2:
############
# Time: 
# Purpose: Trim the data either using BRAT trim function or cut fixed number of bases on each end.
system("date;");
print "------------------------------------------------------------------------\n";
print "Starting Step 2: Trim\n";
print "Before trimming off the low-quality bases, the user can choose to trim the adapter first. The adapter sequence needs to be provided for adapter trimming \n"; 
print "When trimming using BRAT trim function, if the quality score -q value is set up very high (e.g., -q 20), may end up with no reads left \n"; 
print "------------------------------------------------------------------------\n\n";

my $pref = $prefix;



my $retrieve_fixed_trim_fastqOut = "fastX.trimmer.output.fastq";

my $fastx_quality_trim_out = "fastX.quality.trimmer.out";


my $fixedTrim_fastqOut = "fixed.trim.output.fastq";
my $fixedTrim_BRATout = "fixed.brat.compatible.out";



my $trim_input = $input;


	if ($adapter_trim eq "fastx") 
	{
		if ($fastq_format eq "sanger") # Note: this is needed because fastx_clipper only deals with Illumina fastq format
		{
			my $illuminafastq = "illumina.fastq";
			system("perl $code_directory/fq_all2std.ver2.pl sanger2illum $input > $illuminafastq");
			$trim_input = $illuminafastq;
		}
		system("$fastx_dir -a $adapter_sequence -n -l 20 -c -i $trim_input -o fastx.trim.fastq");
		$trim_input = "fastx.trim.fastq";
		# a => ADAPTER string. default is CCTTAAGG (dummy adapter)
		# n => keep sequences with unknown (N) nucleotides. default is to discard such sequences.
		# l => discard sequences shorter than N nucleotides. default is 5.
		# c => Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter)
	}
	if ($adapter_trim eq "cutadapt") 
	{
		system("$cutadapt_dir -a $adapter_sequence $trim_input -m 20 -o cutadapt.trim.fastq");
		$trim_input = "cutadapt.trim.fastq";
		# a => Sequence of an adapter that was ligated to the 3' end.
		# m => Discard trimmed reads that are shorter than LENGTH.
	}

	
	if ($trim_Flag eq "brat") 
	{
   		# This flag means to use BRAT trim function (not "FIX", i.e., not to trim out fixed number of bases)
   		my $trim_quality = "";
   		if ($fastq_format eq "sanger") 
  		{
        	$trim_quality = " -q $brat_trim_qc -L 33 -m $brat_trim_N"; 
   		}
   		elsif ($fastq_format eq "illumina") 
   		{ 
   			$trim_quality = " -q $brat_trim_qc -L 64 -m $brat_trim_N"; 
   		}
       	else 
       	{ 
       		print "The fastq format is neither sanger nor illumina; please double-check the global variables and your fastq data\n";
    	}
   		# q => represents the quality score threshold of bases; bases are trimmed untile a base with a quality score greater or equal than q is encountered
   		# L => specifies the smallest value of the range of base quality scores in ASCII representation. L 33 for sanger and L 64 for Illumina.
   		# m => number of allowable internal Ns

   		system("$brat_trim_dir -s $trim_input -P $pref  $trim_quality");
	}

	if ( $trim_Flag eq "fix" ) 
	{  # This flag means to trim out fixed number of bases set by the users. It is good to use this 
	   # if the users know exactly how many bases they want to trim off and no adaptor sequences 
	   # have been trimmed yet. However, if the users have already trimmed the adapter sequences, 
	   # then the reads length varies, it is NOT good to use this trimming option any more.
  		 system("perl $code_directory/Feb1.trim.fastq.n.bases.pl $trim_input $fixed_trim_5 $fixed_trim_3 $fixedTrim_fastqOut $fixedTrim_BRATout");
	}

  if ($trim_Flag eq "retrieve")
  {
      my $keepBase=$fixed_trim_5+1;
      if ($fastq_format eq "sanger") 
      {
          system("$fastx_trim_dir -Q 33 -f $keepBase -i $trim_input -o $retrieve_fixed_trim_fastqOut");

          system("$fastx_quality_trim_dir -Q 33 -t 20 -l 24 -i $retrieve_fixed_trim_fastqOut -o $fastx_quality_trim_out ");
          system("perl $code_directory/Feb1.trim.fastq.n.bases.pl $fastx_quality_trim_out 0 0 $fixedTrim_fastqOut $fixedTrim_BRATout");

      }

          elsif ($fastq_format eq "illumina") 
      {
          system("$fastx_trim_dir -f $keepBase -i $trim_input -o $retrieve_fixed_trim_fastqOut");

          system("$fastx_quality_trim_dir  -t 20 -l 24 -i $retrieve_fixed_trim_fastqOut -o $fastx_quality_trim_out ");
          system("perl $code_directory/Feb1.trim.fastq.n.bases.pl $fastx_quality_trim_out 0 0 $fixedTrim_fastqOut $fixedTrim_BRATout");

      } 

 
  }

	if ($trim_Flag eq "no") 
	{    # This is the case that the user does not want to do any quality trimming. 
  		 system("perl $code_directory/Feb1.trim.fastq.n.bases.pl $trim_input 0 0 $fixedTrim_fastqOut $fixedTrim_BRATout");
	}


############
# Step 3:
############
# Time: 
# Purpose: brat alignment and acgt count
system("date;");

my $bratIN = $pref . '_reads1.txt';
my $file = "alignment.brat"; # the BRAT alignment output
my $name = "brat.file.name"; # The file contains the path and filename information to obtain the ACGT count


print "------------------------------------------------------------------------\n";
print "Starting Step 3.1: BRAT alignment\n";
print "The users may change the BRAT alignment parameters based on their own data. Currently it is set for single-end reads \n"; 
print "------------------------------------------------------------------------\n\n";
if ($trim_Flag eq "brat") 
{
    print "Now starting to use BRAT to do alignment\n";
    system("$brat_large_dir -r $reference  -s  $bratIN  -S  -m $brat_mismatch  -bs -f $seed -o $file -u -M > brat.log");
}
if (($trim_Flag eq "fix") || ($trim_Flag eq "retrieve")) 
{
    print "Now starting to use BRAT to do alignment\n";
    system("$brat_large_dir -r $reference  -s  $fixedTrim_BRATout  -S  -m $brat_mismatch  -bs -o $file -u -M > brat.log");
}
if ($trim_Flag eq "no") 
{
    print "Now starting to use BRAT to do alignment\n";
    system("$brat_large_dir -r $reference  -s  $fixedTrim_BRATout  -S  -m $brat_mismatch  -bs -o $file -u -M > brat.log");
}
# s => to specify the file with input reads for single-end reads mapping
# S => to specify speed mode for BRAT-large. Space usage with this option doubles, but running time is about three times faster
# bs => to specify bisulfite sequenced reads mapping
# u => is used when a user wishes to output unmapped reads
# M => is used when a user wishes to output ambiguous reads


print "------------------------------------------------------------------------\n";
print "Starting Step 3.2: retrieve trimmed off bases\n\n";
print "------------------------------------------------------------------------\n\n";
if ($trim_Flag eq "retrieve") 
{
    system("perl $code_directory/C.retrieve.pl  $file  $input $reference_directory  $fixed_trim_5  $pref  ");  

    my $retrieved_results = $pref . '.retrieve.alignment'; # output from retrieve step 3.2

    system("echo $retrieved_results >> $name");
    system("date;");
}

else 
{
    system("echo $file >> $name");
    system("date;");
   
}


print "------------------------------------------------------------------------\n";
print "Starting Step 3.3: ACGT-count\n\n";
print "------------------------------------------------------------------------\n\n";
system("$brat_acgt_dir  -r  $reference  -P   $pref  -s $name -B > acgt.log");  

my $sep_IN_forw = $pref . '_forw.txt'; # output from ACGT-count; will be used as input for step 4.1
my $sep_IN_rev = $pref . '_rev.txt';
system("date;");


############
# Step 4.1:
############
system("date;");
print "------------------------------------------------------------------------\n";
print "STEP 4.1: Parsing Data...\n";
print "------------------------------------------------------------------------\n\n";

system("perl $code_directory/data.parser.pl  $sep_IN_forw $output $prefix;");

############
# Step 4.2:
############
# Time: approx <1min
# Purpose: It reads output from step 41. to generate a bisulfate summary table for each chromosome as
#			well as a histogram generated for bisulfate rate for each chromosome.
system("date;");
print "------------------------------------------------------------------------\n";
print "STEP 4.2: Chromosome analysis...\n";
print "------------------------------------------------------------------------\n\n";

system("R CMD BATCH '--args $output $pref ' $code_directory/all.chr.analysis.v3.R");

print "OUTPUT:\n";
print "Chromosome summary table: $output/summary.table.txt\n";
print "Chromosome summary table: $output/CG.coverage.table.txt\n";
print "Bisulfite rate histogram: $output/BS.ps\n";





####################################
# remove all middle step files into middle.step folder
####################################
# `mkdir middle.step.output -p`;
#`mv  *reads1.txt.*    middle.step.output`;
#`mv  *mates*          middle.step.output`;
#`rm  *_rev.txt*`;
#`rm  *brat.file.name`;
#`mv  *_pair1.fastq    middle.step.output`;
#`mv  *base            middle.step.output`;


print "PIPELINE COMPLETE\n";
system("date;");




sub printUsage
{
  print "\n=======================================================================\n"; 
  print "WARNING: 1 or mroe required paramters are missing \n\n";
  print "\nUSAGE: 7 parameters are required \n\n perl BStools.pl -i <FASTQ_input> -t <TARGET_input> -d <MethyQA_directory> ";
  print " -p <prefix> -R <reference_directory> -r <reference_name> [OPTIONS]\n\n"; 
  print " [-i <file>]\t FASTQ input file  \n";
  print " [-t <file>]\t target input file (i.e., a list of target regions specified for analysis) \n";
  print " [-d <dir>]\t path to MethyQA directory (e.g.,  /home/user/downloads/MethyQA/)  \n";
  # print " [-c <string>]\t chromosome number (i.e. chr1, chr2, chr17, chrX, chrY, etc.)\n";
  print " [-p <string>]\t prefix (i.e., the prefix written to the output file names)\n";
  print " [-R <dir>]\t reference directory (i.e., the directory with the genome reference files) \n";
  print " [-r <file>]\t reference name (i.e., the file name of the reference that the user will use) \n";
  print " [-f <string>]\t FASTQ format (i.e., “sanger” or “illumina”) \n";
  print " [-a <string>]\t adapter trimming (i.e., the  user can choose to utilize fastx adapter trimming or cutadapt adapter trimming. Default is no adapter trimming) \n";
  print " \t\t -- no:  no adapter trimming (default) \n";
  print " \t\t -- fastx: fastx adapter trimming \n";
  print " \t\t -- cutadapt: cutadapt adapter trimming. If cutadapt is set, -Y need to be specified in command line \n";
  print " [-A <string>]\t adapter sequence and the default is Illumina adapter sequence AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  \n";
  print " [-T <string>]\t quality trim flag \n";
  print " \t\t -- no:  no quality trimming \n";
  print " \t\t -- brat: brat dynamic trimming (default) \n";
  print " \t\t -- fix: fixed quality trimming \n";
  print " \t\t -- retrieve: fixed quality trimming + retrieve trimmed C sites \n";
  print " [-q <int>]\t for brat trimming, specifies the threshold for quality score";
  print " [-m <int>]\t for brat trimming, specifies the specifies the allowed number of Ns";  
  print " [-N <int>]\t for fixed quality trimming (the users specifies the number of bases to be trimmed at 5' end, default is 5) \n";
  print " [-n <int>]\t for fixed quality trimming (the users specifies the number of bases to be trimmed at 3' end, default is 10) \n";
  print " [-H <int>]\t for brat alignment, specifies the threshold for mismatch";  
  print " [-B <real>]\t cutoff value for high bisulfite rate (([0, 1], default is 0.99))\n";
  print " [-b <real>]\t cutoff value for low bisulfite rate (([0, 1], default is 0.6) \n";
  print " [-L <real>]\t cutoff value for high coverage (([0, 1], default is 0.5) \n";
  print " [-l <real>]\t cutoff value for low coverage (([0, 1], default is 0.1) \n";
  print " [-u <logic>]\t bisulfite flag (it is an option to initiate boxplot of high vs. low bisulfite rates, default is 'TRUE') \n";
  print " [-v <logic>]\t coverage flag (it is an option to initiate boxplot of high vs. low coverage, and default is 'TRUE') \n";
  print " [-Y <string>]\t Pathway to python for running cutadapt(i.e., python, python2.6)\n";
  print " [-Q <string>]\t Pathway to FastQC (i.e., /home/appl/apps/bin/fastqc, default is using the one complied in MethyQA pipeline )\n";
  print " [-M <string>]\t Pathway to BRAT_trim function (i.e., /home/appl/apps/bin/trim.v1.2.4,  default is using the one complied in MethyQA pipeline  )\n";
  print " [-K <string>]\t Pathway to BRAT_large function (i.e., /home/appl/apps/bin/brat-large.v1.2.4, default is using the one complied in MethyQA pipeline  )\n";
  print " [-J <string>]\t Pathway to BRAT_acgt_count function (i.e., /home/appl/apps/bin/brat.v1.2.4, default is using the one complied in MethyQA pipeline  )\n";
  print " [-X <string>]\t Pathway to fastx function (i.e., /home/appl/apps/bin/fastx, default is using the one complied in MethyQA pipeline  )\n";
  print " [-C <string>]\t Pathway to fastx function (i.e., /home/appl/apps/bin/fastx, default is using the one complied in MethyQA pipeline  )\n";
  print  "\n=======================================================================\n\n"; 
}




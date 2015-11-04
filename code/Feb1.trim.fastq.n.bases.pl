# File name:  Feb1.trim.fastq.n.bases.pl
# Purpose:    trim reads for last n bases at 3' end and m bases at 5'end 
# Example:    Feb1.trim.fastq.n.bases.pl $ARGV[0]  5   10  fixedTrim_fastqOut  fixedTrim_BRATout
# Time:       It takes about 3 minutes to trim 26 million of reads. 

use strict; use warnings;

print "\n--This perl code to trim reads for first m bases at 5' end and last n bases at 3' end . \n--There are  five arguments:\n"; 
print "--(1) the fastq file, (2) number of bases to trim at the 5' end (i.e. the left side), (3) number of bases to trim at the 3' end (i.e. the right side) and (4) output file name. (5) second output file\n--The usage is: perl  Feb1.trim.fastq.n.bases.pl  Feb1.5reads  left.mbase  right.nbase  Feb.5read.trim\n\n"; 

open IN, "$ARGV[0]" || die "Cannot find file$!";
open OUT, ">$ARGV[3]" || die "Cannot find file$!";
open OUT1, ">$ARGV[4]" || die "Cannot find file$!"; 

# Define gloable variables:

my $read; my $quality;  my $third_line;  my $new_read; 
my $new_quality;  my $read_length; my $count = 0; 

# First we check the length of a read by checking the first read

while(<IN>) 
{
   while($_ =<IN>) # "$_" contains each line of the file. 
   { 
     chomp; 
     $read_length = length($_);
   }
}
seek(IN, 0, 0); # to make sure next we will start from the begining to the end

if ($read_length <= ($ARGV[1]+$ARGV[2]))
{   print ("read length is $read_length. It is less or equal to $ARGV[1]+$ARGV[2]. You are cutting off too many bases. \n")  }
else 
{ 
  # seek(IN, 0, 0); # to make sure next we will start from the begining to the end
  # Note, I should not have the above seek function inside of "else". It gives wrong results. 
  # I should move it out side of this "if ... else" condition. 

  while(<IN>)
  { 
     if (/^@/) 
     { print OUT;                  # print out the header line

       $read = <IN>;  chomp($read); 
       $new_read=substr($read, $ARGV[1], $read_length - $ARGV[1] - $ARGV[2]); 
       print OUT "$new_read\n";    # Print out the trimmed read
	   print OUT1 "$new_read\t$ARGV[1]\t$ARGV[2]\n"; # second output: sequence \t bases trimmed at 5' \t bases trimmed at 3'
       
       $third_line=<IN>; chomp($third_line); 
       print OUT "$third_line\n";  # print out the header line again

       $quality=<IN>; chomp($quality); 
       $new_quality = substr($quality, $ARGV[1], $read_length  - $ARGV[1] - $ARGV[2]); 
       print OUT "$new_quality\n"; # Print out the trimmed quality score
    }
     
     # Print time: 
     $count = $count + 1;  # print "count is: $count \n"; 
     if ( $count%1000000 ==0 ) 
     {  my $now = localtime time;
        print ("-- When finishing $count reads, time is: $now. \n") ; 
     }
  }
}
      
print "-------------------------- Finished trimming \n" ; 

close IN;
close OUT;
close OUT1;


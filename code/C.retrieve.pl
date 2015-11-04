#!/usr/bin/perl -w

use warnings;
use strict;


#file usage, checks for valid number of input parameters
my $usage = qq(Usage:  perl C.retrive.pl </path/to/alignment/file>  <path/to/raw/sequencing/file> <path/to/reference/directory> <retrieved base> <output name>);

die($usage) if (@ARGV < 5);

#create the file handlers and open

open(OUTPUT, ">", "$ARGV[4].retrieve.alignment") or die $!;	 #the date will be attached to the prefix of the output

my $Ref_folder = $ARGV[2];
my $retrieved_base = $ARGV[3];

#global headers
my $alignment; my @alignmentSplit;
my $index; my $read; my $read_length; my $chr; my $strand; my $pos; my $mismatch; my $rawPos;
my $ref_file; my $start_pos; 
my $retrive_seq; my $retrive_1base;
my $new_pos; 
my $count=0;
my %hash;  my $raw_read; my $third_line; my $quality; my $matched_raw_read;


open(Sequence, $ARGV[1]) || die("the path is incorrect") or die $!; #$ARGV[1]
  while(<Sequence>)
  { 
     
     if (/^@/) 
     { 
       $raw_read = <Sequence>;  chomp($raw_read); 
       $third_line=<Sequence>; chomp($third_line); 
       $quality=<Sequence>; chomp($quality); 

    }

       $count = $count + 1;  # print "count is: $count \n"; 
       $hash{$count} = $raw_read;
 }


open(Alignment, $ARGV[0]) || die("the path is incorrect") or die $!; #$ARGV[1]

while($alignment = <Alignment>){ #while = for each aligned read, do the following:
	chomp $alignment;
	@alignmentSplit = split("\t", $alignment);
	$index = $alignmentSplit[0];
	$read = $alignmentSplit[1];
	$read_length =length($read);
	$chr = $alignmentSplit[2];
	$strand = $alignmentSplit[3];
	$pos = $alignmentSplit[4];	
	$mismatch = $alignmentSplit[5];	
	$rawPos = $alignmentSplit[6];

	$matched_raw_read = $hash{$index+1};


    $ref_file = $Ref_folder.'/'.$chr.'.fa';

    if ($strand eq "+"){
         $start_pos=$pos+1-5;
         $retrive_seq = obtain_subroutine($ref_file, $start_pos, $retrieved_base);
         $retrive_1base = substr($retrive_seq,0,1);
         if ($retrive_1base eq "C" ||  $retrive_1base eq "c"){
         	   $new_pos = $pos-$retrieved_base;
         	   print (OUTPUT "$index\t$matched_raw_read\t$chr\t$strand\t$new_pos\t$mismatch\t$rawPos\n");
         }
    }
    else {
    	 $start_pos=$pos+1+$read_length;
         $retrive_seq = reverse(obtain_subroutine($ref_file, $start_pos, $retrieved_base));
    	 $retrive_seq =~ tr/ACGTacgt/TGCAtgca/;
         $retrive_1base = substr($retrive_seq,0,1);
         if ($retrive_1base eq "C" ||  $retrive_1base eq "c"){
              $new_pos = $pos;
         	   print (OUTPUT "$index\t$matched_raw_read\t$chr\t$strand\t$new_pos\t$mismatch\t$rawPos\n");
         }
    }

}

close Alignment;
close Sequence;
close OUTPUT;




sub obtain_subroutine
{
my ($ref, $start, $retrive_length) = @_; ## this holds the file name passed when the subroutine is called
open (CHR, "<$ref") or die "error reading file $ref";

 my $firstline=<CHR>;
       my $first_line_length=length($firstline);
        chomp($firstline);
       my $line_length=length($firstline);
       my @chr_name=split(">",$firstline); ##This is done to obtain the chrname from the input FASTA file

       #####################################################################################
        my $a=int($start/50);
        my $b=$start%50;

		if($a>0 && $b>0) ### This is done to ignore the /n characters residing inside the FASTA files
                { #print("here1");
                $start+=$a;
                }
       if($a>0 && $b==0) ### This is done to ignore the /n characters residing inside the FASTA files
       {
                #print("here2");
                $start=$start+$a-1;
                }

		my $origin_seq;
        for (my $i=$start;$i<=$start+$retrive_length-1;$i++)       ## All of the lines in this section are used to obtain the original sequence from the FASTA
        {                                                       ## files from the start and end postion mention by the user.
                seek(CHR,$i+$first_line_length-1,0);
                my $mk=getc(CHR);
                if($mk eq "\n")
                {
                        $retrive_length++;
                }
                 else
                {
                        $origin_seq=$origin_seq.$mk;
                }

        }
        return $origin_seq;### The origin_seq variable holds the original sequence from the FASTA file
        close CHR;

}

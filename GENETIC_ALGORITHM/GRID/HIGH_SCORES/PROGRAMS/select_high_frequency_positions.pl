#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $inputFilename="";
   my $thresholdFrequency=0;
   my @con=();
   my $line="";
   my $interfacePosition="";
   my $frequency=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----

# Main code of the program starts here.

if($#ARGV < 0)
{
   print "\nUsage: $0 <Name of input file> <Threshold frequency>";
   print "\n\n";
   exit(0);
}

$inputFilename = $ARGV[0];
$thresholdFrequency = $ARGV[1];

if(! -e $inputFilename)
{
   print "\nFile \"$inputFilename\" does not exist. Aborting program";
   exit(0);
}

open(hd,$inputFilename);
@con=<hd>;
close(hd);

chomp(@con);

foreach $line (@con)
{
   $line=~s/^ +//g;
   $line=~s/\t+/ /g;
   $line=~s/\s+/ /g;

   ($interfacePosition,$frequency)=split(/\s/,$line);

   if($frequency >= $thresholdFrequency)
   {
      $interfacePosition=sprintf("%8s",$interfacePosition);
      $frequency=sprintf("%10d",$frequency);
      print $interfacePosition,$frequency,"\n";
   }
}

# End of program.

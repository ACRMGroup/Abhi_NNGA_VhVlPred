#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $inputFilename="";
   my @con=();
   my $line="";
   my $ignore="";
   my $predicted=0;
   my $actual=0;
   my $error=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----

# Main code of the program starts here.

if($#ARGV < 0)
{
   print "\nUsage: $0 <Name of file with the predicted and actual output values>\n\n";
   exit(0);
}

$inputFilename=$ARGV[0];

if(! -e $inputFilename)
{
   print "\nFile \"$inputFilename\" does not exist.\nAborting program\n\n";
   exit(0);
}

open(hd,$inputFilename);
@con=<hd>;
close(hd);

@con=grep(/Predicted/,@con);

chomp(@con);

foreach $line (@con)
{
   chomp($line);

   # (Predicted,Actual): 0.5268 0.5633

   ($ignore,$line)=split(/:/,$line);
   $line=~s/^\s//;

   ($predicted,$actual)=split(/ /,$line);

   $error=$actual-$predicted;

   print $error,"\n";
}

# End of program.

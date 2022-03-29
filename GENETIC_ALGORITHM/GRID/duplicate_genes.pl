#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $inputFilename="";
   my $geneNumber=0;
   my $outputFilename="";
   my $numberOfGenes=0;
   my $generation=0;
   my $appendMode="Y";

   my @con="";
   my $startString="";
   my $i=0;

   my @currentGenes=();
   my $numberOfGenesAlreadyPresent=0;
   my $geneNumberStart=0;
   my $geneNumberEnd=0;

   my $ignore="";
   my $gene="";
   my $score=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# Main code of the program starts here.

if($#ARGV < 5)
{
   print "\nUsage: $0 <Arguments>\n";
   print "\nArguments are:\n";
   print "\n1. Name of file with the genes";
   print "\n2. Gene number";
   print "\n3. Name of output file";
   print "\n4. Number of genes to be duplicated";
   print "\n5. Generation";
   print "\n6. Append genes to the output file - Y or N";
   print "\n\n";
   exit(0);
}

$inputFilename=$ARGV[0];
$geneNumber=$ARGV[1];
$outputFilename=$ARGV[2];
$numberOfGenes=$ARGV[3];
$generation=$ARGV[4];
$appendMode=$ARGV[5];

open(hd,$inputFilename);
@con=<hd>;
close(hd);

$startString="GENE #".$geneNumber;

$i=0;

while($con[$i] !~ $startString)
{
   $i++;
}

# Line resembles GENE #1: 0001000010100010001001110000000.
# Extract the gene.

($ignore,$gene)=split(/:/,$con[$i]);
$gene=~s/\s//g;

# Pull out the score for the gene.

$i++;

($ignore,$score)=split(/:/,$con[$i]);
$score=~s/\s//g;

# Write the gene into the output file.

if($appendMode ne "N")
{
   # Read the existing content to ascertain how many genes are already present in the file.

   open(hd,$outputFilename);
   @con=<hd>;
   close(hd);

   @currentGenes=grep(/GENE /,@con);
   $numberOfGenesAlreadyPresent=$#currentGenes+1;

   $geneNumberStart=$numberOfGenesAlreadyPresent + 1;
   $geneNumberEnd=$geneNumberStart + $numberOfGenes;

   open(whd,">>$outputFilename");
}
else
{
   $geneNumberStart=1;
   $geneNumberEnd=$numberOfGenes;

   open(whd,">$outputFilename");
}

for($i=$geneNumberStart;$i < $geneNumberEnd;$i++)
{
   print whd "GENE #",$i,": $gene\nSCORE: $score\nGENERATION: $generation\n";
   print whd "-------------------------------------\n";
}

close(whd);

# End of program.

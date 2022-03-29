#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $inputFilename="";
   my @con=();
   my @existingScores=();
   my $numberOfGenesAlreadyPresent=0;
   my $line="";
   my $gene="";
   my $score=0;
   my @parts=();
   my %hash=();
   my @scores=();
   my @sortedScores=();
   my $outputFilename="";
   my $i=0;
   my $geneNumberStart=0;
   my $geneNumberEnd=0;
   my $j=0;
   my $generation=0;
   my $N=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# Main code of the program starts here.

if($#ARGV < 3)
{
   print "\nUsage: $0 <Input file> <Number of genes to write> <Output file> <Generation>";
   print "\n\n";
   exit(0);
}

$inputFilename=$ARGV[0];
$N=$ARGV[1];
$outputFilename=$ARGV[2];
$generation=$ARGV[3];

# Read existing scores into an array.

open(hd,$outputFilename);
@con=<hd>;
close(hd);

@existingScores=grep(/SCORE/,@con);

$numberOfGenesAlreadyPresent=$#existingScores+1;
$geneNumberStart=$numberOfGenesAlreadyPresent+1;
$geneNumberEnd=$geneNumberStart + $N;

# Read contents of current file.

open(hd,$inputFilename);
@con=<hd>;
close(hd);

foreach $line (@con)
{
   chomp($line);

   if($line =~ /GENE /)
   {
      $gene=$line;
      $gene=~s/GENE.*: //;
   }
   elsif($line =~ /SCORE/)
   {
      $score=$line;
      $score=~s/SCORE: //;

      @parts=grep(/$score/,@scores);

      if($#parts != -1)
      {
	 next;
      }

      @parts=grep(/$score/,@existingScores);

      if($#parts != -1)
      {
	 next;
      }

      $hash{$score}=$gene;
      push(@scores,$score);
   }
}

# Sort the scores.

@sortedScores = sort {$a <=> $b} @scores;

# Write the N best scores into the output file.

open(whd,">>$outputFilename");

for($i=$#sortedScores,$j=$geneNumberStart;($j < $geneNumberEnd) && ($i >= 0);$i--,$j++)
{
   print whd "GENE #$j: ",$hash{$sortedScores[$i]},"\n","SCORE: $sortedScores[$i]\n";
   print whd "GENERATION: $generation\n";
   print whd "---------------------------\n";
}

close(whd);

# End of program.

#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my @con=();
   my $numberOfRedundantSequences=0;

   my @generations=0;
   my @redundant=();
   my @bestScores=();
   my @averageScores=();
   my @standardDeviations=();

   my $i=0;
   my $printHeader=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# --------------- SUB - ROUTINES SECTION ---------------

sub get_redundancy
{
   my $numberOfLines=$#con+1;
   my $i=0;
   my $line="";
   my @allUniqueGenes=();
   my $numberOfGenes=0;
   my @parts=();

   for($i=0;$i < $numberOfLines;$i++)
   {
      $line=$con[$i];

      if($line =~ /GENE #/)
      {
	 $line=~s/GENE #.*://;
	 $line=~s/\s//g;

	 @parts=grep(/$line/,@allUniqueGenes);	

	 if($#parts == -1)
	 {
	    push(@allUniqueGenes,$line);
	 }

	 $numberOfGenes++;
      }
      elsif( ($line =~ /GENERATION #/) && ($#allUniqueGenes != -1) )
      {
	 $numberOfRedundantSequences=$numberOfGenes - ($#allUniqueGenes + 1);
	 push(@redundant,$numberOfRedundantSequences);

	 $numberOfRedundantSequences=0;
	 $numberOfGenes=0;
	 @allUniqueGenes=();
      }
   }

   # Calculate statistic for the last generation.

   $numberOfRedundantSequences=$numberOfGenes - ($#allUniqueGenes + 1);
   push(@redundant,$numberOfRedundantSequences);

} # End of function "get_redundancy".


sub get_score_stats
{
   my $numberOfLines=$#con+1;
   my $i=0;
   my $line="";
   my $numberOfGenes=0;

   my @allScores=();
   my $bestScore=0;
   my $averageScore=0;
   my $standardDeviation=0;

   for($i=0;$i < $numberOfLines;$i++)
   {
      $line=$con[$i];

      if($line =~ /SCORE/)
      {
	 $line=~s/SCORE: //;

	 push(@allScores,$line);
      }
      elsif( ($line =~ /GENERATION #/) && ($#allScores != -1) )
      {
	 # Calculate all the numbers required.

	 ($bestScore,$averageScore,$standardDeviation)=&calculate_stats(\@allScores);

	 push(@bestScores,$bestScore);
	 push(@averageScores,$averageScore);
	 push(@standardDeviations,$standardDeviation);

	 @allScores=();
      }
   }

   # Calculate numbers for the last generation.

   ($bestScore,$averageScore,$standardDeviation)=&calculate_stats(\@allScores);

   push(@bestScores,$bestScore);
   push(@averageScores,$averageScore);
   push(@standardDeviations,$standardDeviation);

} # End of function "get_score_stats".


sub calculate_stats
{
   my $scoresArgument=$_[0];
   my @scores=@$scoresArgument;
   my $i=0;
   my $bestScore=0;
   my $averageScore=0;
   my $standardDeviation=0;
   my $variance=0;
   my $total=0;

   my $numberOfScores=($#scores + 1);

   for($i=0;$i <= $numberOfScores;$i++)
   {
      $total+=$scores[$i];

      if($bestScore < $scores[$i])
      {
	 $bestScore=$scores[$i];
      }
   }

   $averageScore=$total/$numberOfScores;

   $total=0;

   for($i=0;$i < $numberOfScores;$i++)
   {
      $total = $total + ( ($scores[$i] - $averageScore) * ($scores[$i] - $averageScore) );
   }

   $variance=$total/$numberOfScores;
   $standardDeviation=sqrt($variance);

   return $bestScore,$averageScore,$standardDeviation;

} # End of function "calculate_stats".


# ----------- END OF SUB - ROUTINES SECTION ------------

# Main code of the program starts here.

if($#ARGV < 0)
{
   print "\nUsage: $0 <Input file> <-h (Optional - to print header before table)>\n\n";
   exit(0);
}

open(hd,$ARGV[0]);
@con=<hd>;
close(hd);

if($ARGV[1] eq "-h")
{
   $printHeader=1;
}
else
{
   $printHeader=0;
}

# Get generations.

@generations=grep(/GENERATION #/,@con);
chomp(@generations);

for($i=0;$i <= $#generations;$i++)
{
   $generations[$i]=~s/GENERATION #//;
}

# Find number of redundant sequences.

&get_redundancy();

# Get the overall score statistics.

&get_score_stats();

# Remove newline characters from all arrays.

chomp(@generations);
chomp(@redundant);
chomp(@bestScores);
chomp(@averageScores);
chomp(@standardDeviations);

# Print the values.

if($printHeader)
{
   print "\n\n";
   print "Generation     Redundant\tBest-Score\tAverage-Score\t\tSigma\n";
   print "               sequences\n\n";
}

for($i=0;$i <= $#generations;$i++)
{
   printf "%5d",$generations[$i];
   printf "%15d",$redundant[$i];
   printf "%20f",$bestScores[$i];
   printf "%20f",$averageScores[$i];
   printf "%20f",$standardDeviations[$i];
   print "\n";
}

print "\n";

# End of program.

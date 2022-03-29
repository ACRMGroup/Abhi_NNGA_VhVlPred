#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $genesFilename="";
   my $interfacePositionsListFilename="";
   my @con=();
   my @genes=();
   my $line="";
   my @interfacePositionsList=();
   my $gene="";
   my $geneNumber=0;
   my $geneString="";
   my @zeroesAndOnes=();
   my $i=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# Main code of the program starts here.

if($#ARGV < 1)
{
   print "\nUsage: $0 <File with genes> <File with interface positions list>";
   print "\n\n";
}

$genesFilename=$ARGV[0];
$interfacePositionsListFilename=$ARGV[1];

# Read all the genes into an array.

open(hd,$genesFilename);
@con=<hd>;
close(hd);

@genes=grep(/GENE #/,@con); # GENE #1: 00000001000000000000100110010000000000000001000000010010

chomp(@genes);

# Read the interface positions list.

open(hd,$interfacePositionsListFilename);
@con=<hd>;
close(hd);

foreach $line (@con)
{
   if($line =~ /^L/)
   {
      push(@interfacePositionsList,$line);
   }
   elsif($line =~ /^H/)
   {
      push(@interfacePositionsList,$line);
   }
}

chomp(@interfacePositionsList);

# Now, correlate the interface positions with the genes.

foreach $gene (@genes)
{
   # Isolate the gene string from "GENE #1: 00000001000000000000100110010000000000000001000000010010"

   ($geneNumber,$geneString)=split(/: /,$gene);

   # Split the gene string into 0s and 1s.

   @zeroesAndOnes=split(//,$geneString);

   # Now, write out the list of interface positions.

   print ">$geneNumber\n";

   for($i=0;$i <= $#zeroesAndOnes;$i++)
   {
      if($zeroesAndOnes[$i] == 1)
      {
	 print $interfacePositionsList[$i],"\n";
      }
   }
}

# End of program.

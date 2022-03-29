#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $limitedInterfacePositionsFilename="";
   my $allInterfacePositionsFilename="";

   my @con=();

   my @lightChainPositions=();
   my @heavyChainPositions=();
   my @limitedInterfacePositions=();
   my @allInterfacePositions=();

   my $i=0;
   my $interfacePosition="";
   my @geneString=();
   my @parts=();

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----

if($#ARGV < 1)
{
   print "\nUsage: $0 <Arguments>";
   print "\n\n";
   print "Arguments are:\n\n";
   print "1. File with the limited set of interface positions\n";
   print "2. File with all the interface positions";
   print "\n\n";
   exit(0);
}

$limitedInterfacePositionsFilename = $ARGV[0];
$allInterfacePositionsFilename = $ARGV[1];

# Read the file with the limited set of interface positions.

open(hd,$limitedInterfacePositionsFilename);
@con=<hd>;
close(hd);

@lightChainPositions=grep(/^L[0-9]/,@con);
@heavyChainPositions=grep(/^H[0-9]/,@con);

@limitedInterfacePositions=@lightChainPositions;
push(@limitedInterfacePositions,@heavyChainPositions);

undef(@con);
undef(@lightChainPositions);
undef(@heavyChainPositions);

# Read file with all interface positions.

open(hd,$allInterfacePositionsFilename);
@con=<hd>;
close(hd);

@lightChainPositions=grep(/^L[0-9]/,@con);
@heavyChainPositions=grep(/^H[0-9]/,@con);

@allInterfacePositions=@lightChainPositions;
push(@allInterfacePositions,@heavyChainPositions);

# Construct the gene.

@geneString="";

print "GENE #1: ";

for($i=0;$i <= $#allInterfacePositions; $i++)
{
   $interfacePosition=$allInterfacePositions[$i];

   @parts=grep(/($interfacePosition)$/,@limitedInterfacePositions);

   if($#parts != -1)
   {
      push(@geneString,'1');
   }
   else
   {
      push(@geneString,'0');
   }

   print $geneString[$i];
}

print "\nSCORE: 0\nGENERATION: 1\n";
print "-------------------------------------\n";

# End of program.

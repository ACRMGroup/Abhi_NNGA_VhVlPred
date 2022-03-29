#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $untrainedNetworkPath="";
   my $intermediateFilename="";
   my $numberOfTrainingPatterns=400;
   my $targetDirectory="";
   my $untrainedNetworkType="";

   my @con=();
   my @pdbCodes=();
   my $numberOfPDB=-1;

   my $boundaryStart=0;
   my $boundaryEnd=0;
   my $trainingSetRange="";
   my $outputFilename="";

   my $scriptFilename="";

   my $HOME="";
   my $command="";
   my $snnsConfigFilename="";

   my @interfacePositionsMethods=("frequency_average","frequency_max_min","max_min","average");
   my $interfacePositionsMethodsIndex=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# --------------- SUB - ROUTINES SECTION ---------------

sub Usage
{
   print "\n./$0 <Arguments>\n";
   print "\nArguments are:\n";
   print "\n1. Path of untrained network";
   print "\n2. Path of intermediate file";
   print "\n3. Target directory for script";
   print "\n4. Number of training set patterns (Default: 400) - Optional parameter";
   print "\n5. Path of file to configure parameters for SNNS - Optional parameter";
}

# ----------- END OF SUB - ROUTINES SECTION ------------


# Main code of the program starts here.

# Format of the program is as follows:
#
# ./final.exe <Arguments>
#
# -net UNTRAINED_NEURAL_NETWORKS/20i_20h.net
# -trb 1-400
# -int INTERMEDIATE_FILES/frequency_average_intermediate.dat
# -out 1-400.out

if($#ARGV < 2)
{
   &Usage();
   print "\n\n";
   exit(0);
}

$untrainedNetworkPath=$ARGV[0];
$intermediateFilename=$ARGV[1];
$targetDirectory=$ARGV[2];

if($#ARGV >= 3)
{
   $numberOfTrainingPatterns=$ARGV[3];
}

if(-e $ARGV[4])
{
   $snnsConfigFilename=$ARGV[4];
}

open(hd,$intermediateFilename);
@con=<hd>;
close(hd);

$HOME=$ENV{"HOME"};

@pdbCodes=grep(/PDB Code/,@con);
$numberOfPDB=$#pdbCodes+1;

@con=split(/\//,$untrainedNetworkPath);
$con[$#con]=~s/\.net//;
$untrainedNetworkType=$con[$#con];

@con=();

$boundaryStart=1;
$boundaryEnd=$numberOfTrainingPatterns;

while($interfacePositionsMethodsIndex <= $#interfacePositionsMethods)
{
   if($intermediateFilename =~ /$interfacePositionsMethods[$interfacePositionsMethodsIndex]/)
   {
      last;
   }

   $interfacePositionsMethodsIndex++;
}

if($interfacePositionsMethodsIndex <= $#interfacePositionsMethods)
{
   $scriptFilename=$targetDirectory."/script_".$interfacePositionsMethods[$interfacePositionsMethodsIndex];
   $scriptFilename.="_".$untrainedNetworkType."_".$numberOfTrainingPatterns.".sh";
}
else
{
   $scriptFilename=$targetDirectory."/script_".$untrainedNetworkType."_".$numberOfTrainingPatterns.".sh";
}

$outputFilename=$numberOfTrainingPatterns."_".$untrainedNetworkType.".out";

if($untrainedNetworkPath !~ /^\//)
{
   $untrainedNetworkPath="$HOME/SNNS/$untrainedNetworkPath";
}

if($intermediateFilename !~ /^\//)
{
   $intermediateFilename="$HOME/SNNS/$intermediateFilename";
}

open(whd,">$scriptFilename");

print whd "#! /usr/bin/sh\n\n";
print whd "# Number of training patterns: $numberOfTrainingPatterns\n\n";
print whd "rm -f $outputFilename\n\n";

if(-e $snnsConfigFilename)
{
   print whd <<EOT;

echo >> $outputFilename
echo "Config file: $snnsConfigFilename" >> $outputFilename
echo >> $outputFilename

EOT
}

while($boundaryStart != $numberOfTrainingPatterns)
{
   $trainingSetRange=$boundaryStart."-".$boundaryEnd;

   $command="$HOME/SNNS/final.exe -net $untrainedNetworkPath";
   $command.=" -trb $trainingSetRange";
   $command.=" -int $intermediateFilename";
   $command.=" -out $outputFilename";

   if(-e $snnsConfigFilename)
   {
      $command.=" -scfg $snnsConfigFilename";
   }

   print whd $command,"\n\n";

   $boundaryStart++;
   $boundaryEnd++;

   if($boundaryEnd > $numberOfPDB)
   {
      $boundaryEnd=1;
   }
}

close(whd);

# End of program.

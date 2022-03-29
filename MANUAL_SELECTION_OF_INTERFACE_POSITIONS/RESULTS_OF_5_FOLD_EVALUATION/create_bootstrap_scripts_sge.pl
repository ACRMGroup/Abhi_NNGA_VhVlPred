#! /acrm/usr/local/bin/perl

use strict 'vars';
use Cwd;

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $numberOfPDB=0;
   my $genesFilename="";
   my $geneNumber="";
   my $numberOfPatternsPerRun=0;

   my $abInterface="";
   my $programPath="";
   my $pdbListFilename="";
   my $interfaceAnglesFilename="";
   my @pdbList=();
   my $numberOfPatterns=0;
   my $interfacePositionsListFilename="";

   my $boundaryStart=0;
   my $boundaryEnd=0;
   my $trainingSetRange="";
   my $command="";

   my $scriptFilename="";
   my $outputFilename="";
   my $currentWorkingDirectory="";

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# --------------- SUB - ROUTINES SECTION ---------------

sub Usage
{
   print "\n$0 <Arguments>\n";
   print "\nArguments are:\n";
   print "\n1. Path of file containing the gene(s)";
   print "\n2. Number of gene in the file containing the gene(s)";
   print "\n3. Number of training set patterns per run (Default: 400) - Optional parameter";
}


sub write_header
{
   print whd <<EOT;
#!/bin/sh
#
#
#___INFO__MARK_BEGIN__
##########################################################################
#
#  The Contents of this file are made available subject to the terms of
#  the Sun Industry Standards Source License Version 1.2
#
#  Sun Microsystems Inc., March, 2001
#
#
#  Sun Industry Standards Source License Version 1.2
#  =================================================
#  The contents of this file are subject to the Sun Industry Standards
#  Source License Version 1.2 (the "License"); You may not use this file
#  except in compliance with the License. You may obtain a copy of the
#  License at http://gridengine.sunsource.net/Gridengine_SISSL_license.html
#
#  Software provided under this License is provided on an "AS IS" basis,
#  WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING,
#  WITHOUT LIMITATION, WARRANTIES THAT THE SOFTWARE IS FREE OF DEFECTS,
#  MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE, OR NON-INFRINGING.
#  See the License for the specific provisions governing your rights and
#  obligations concerning the Software.
#
#  The Initial Developer of the Original Code is: Sun Microsystems, Inc.
#
#  Copyright: 2001 by Sun Microsystems, Inc.
#
#  All Rights Reserved.
#
##########################################################################
#___INFO__MARK_END__

# This is a simple example of a SGE batch script

# request Bourne shell as shell for job

EOT

print whd "#\$ -S /bin/sh\n\n";

print whd <<EOT;

export LD_LIBRARY_PATH="/usr/local/lib:/acrm/gridengine/lib/grass6.0.0-i686-pc-linux-gnu-14_03_2005/lib/:/home/bsm2/abhi/lib:/home/bsm/martin/lib"

export TEMPORARY_SNNS_PATH="/acrm/home/abhi/tmp/SNNS/"
export ABINTERFACE="/home/bsm2/abhi/VH-VL-INTERFACE/STRUCTURAL_ANALYSIS/"
export ABHIDATADIR="/home/bsm2/abhi/DATA/"
export PATH="/acrm/gridengine/bin/glinux:/usr/kerberos/bin:/acrm/usr/local/bin:/usr/bin:/usr/local/bin:/bin:/usr/X11R6/bin:/home/bsm2/abhi/BASIC_UTILITIES:/home/bsm/martin/bin:/home/bsm2/abhi/VH-VL-INTERFACE/KABAT_NUMBERING/"

EOT

} # End of function "write_header".

# ----------- END OF SUB - ROUTINES SECTION ------------


# Main code of the program starts here.

# Format of the program is as follows:
#
# ./final.exe <Arguments>
#
# -in GENETIC_ALGORITHM/GRID/HIGH_SCORES/RUN2/gene_0.6785.txt
# -num 528
# -int GENETIC_ALGORITHM/GRID/all_interface_positions.dat
# -pdb $ABINTERFACE/Fv_Fab_pdb.lst
# -tor $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt
# -trb 1-480
# -out OUTPUTS/1-480.out

if($#ARGV < 2)
{
   &Usage();
   print "\n\n";
   exit(0);
}

$currentWorkingDirectory=getcwd();

$genesFilename=$ARGV[0];
$geneNumber=$ARGV[1];
$numberOfPatternsPerRun=$ARGV[2];

if($numberOfPatternsPerRun == 0)
{
   $numberOfPatternsPerRun=400;
}

$abInterface=$ENV{"ABINTERFACE"};

$programPath=$abInterface."/final.exe";
$interfacePositionsListFilename=$abInterface."/all_interface_positions.dat";
$pdbListFilename=$abInterface."/Fv_Fab_pdb.lst";
$interfaceAnglesFilename=$abInterface."/INTERFACE_ANGLE/torsion_angles.txt";

open(hd,$pdbListFilename);
@pdbList=<hd>;
close(hd);

$numberOfPatterns=$#pdbList + 1;

$boundaryStart=1;

while($boundaryStart != $numberOfPatterns)
{
   $boundaryEnd=$boundaryStart + $numberOfPatternsPerRun - 1;

   if($boundaryEnd > $numberOfPatterns)
   {
      $boundaryEnd=$boundaryEnd - $numberOfPatterns;
   }

   $trainingSetRange=$boundaryStart."-".$boundaryEnd;

   $outputFilename=$currentWorkingDirectory."/OUTPUTS/$trainingSetRange.out";
   $scriptFilename=$currentWorkingDirectory."/SCRIPTS/s"."$trainingSetRange.csh";

   open(whd,">$scriptFilename");

   &write_header();

   # ./final.exe <Arguments>
   #
   # -in GENETIC_ALGORITHM/GRID/HIGH_SCORES/RUN2/gene_0.6785.txt
   # -num 528
   # -int GENETIC_ALGORITHM/GRID/all_interface_positions.dat
   # -pdb $ABINTERFACE/Fv_Fab_pdb.lst
   # -tor $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt
   # -trb 1-480
   # -out OUTPUTS/1-480.out

   $command=$programPath." -in $genesFilename";
   $command.=" -num $geneNumber";
   $command.=" -int $interfacePositionsListFilename";
   $command.=" -pdb $pdbListFilename";
   $command.=" -tor $interfaceAnglesFilename";
   $command.=" -trb $trainingSetRange";
   $command.=" -out $outputFilename";

   print whd $command,"\n\n";

   $boundaryStart++;
   $boundaryEnd++;

   if($boundaryEnd > $numberOfPDB)
   {
      $boundaryEnd=1;
   }

   close(whd);

   chmod 0755,$scriptFilename;
}

# End of program.

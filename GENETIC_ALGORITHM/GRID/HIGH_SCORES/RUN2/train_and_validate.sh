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
#$ -S /bin/sh

# Parameters to the program:
#
# ./train_and_validate.sh $inputFilename $genePopulation $geneNumber $outputDirectory $generation $geneType $mutationRate

if [ $# -lt 7 ]
then
   echo
   echo "Usage: $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. Input filename"
   echo "2. Gene population"
   echo "3. Gene number"
   echo "4. Output directory"
   echo "5. Generation"
   echo "6. Gene type - P or C"
   echo "7. Mutation Rate"
   echo "8. Number of bins"
   echo
   exit 0
fi

inputFilename=$1
genePopulation=$2
geneNumber=$3
outputDirectory=$4
generation=$5
geneType=$6
mutationRate=$7
numberOfBins=$8

if [ "$numberOfBins" != "" ]
then
   export NUMBER_OF_BINS=$numberOfBins
fi

export LD_LIBRARY_PATH="/usr/local/lib:/acrm/gridengine/lib/grass6.0.0-i686-pc-linux-gnu-14_03_2005/lib/:/home/bsm2/abhi/lib:/home/bsm/martin/lib"

export GENE_POPULATION=$genePopulation
export TEMPORARY_SNNS_PATH="/acrm/home/abhi/tmp/SNNS/"
export ABINTERFACE="/home/bsm2/abhi/VH-VL-INTERFACE/STRUCTURAL_ANALYSIS/"
export ABHIDATADIR="/home/bsm2/abhi/DATA/"
export PATH="/acrm/gridengine/bin/glinux:/usr/kerberos/bin:/acrm/usr/local/bin:/usr/bin:/usr/local/bin:/bin:/usr/X11R6/bin:/home/bsm2/abhi/BASIC_UTILITIES:/home/bsm/martin/bin:/home/bsm2/abhi/VH-VL-INTERFACE/KABAT_NUMBERING/"
export MUTATION_RATE=$mutationRate
export POSITIONS_TO_SWAP=5
export HOSTNAME=`hostname`

echo "MUTATION_RATE: $MUTATION_RATE"

# Execute the program "train_and_validate.exe"

if [ "$numberOfBins" != "" ]
then

   ./train_and_validate.exe -in $inputFilename -num $geneNumber -int interface_positions_best_gene.txt -pdb $ABINTERFACE/Fv_Fab_pdb.lst -tor $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt -type $geneType -odir $outputDirectory -gen $generation -bins

else
 
   ./train_and_validate.exe -in $inputFilename -num $geneNumber -int interface_positions_best_gene.txt -pdb $ABINTERFACE/Fv_Fab_pdb.lst -tor $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt -type $geneType -odir $outputDirectory -gen $generation

fi

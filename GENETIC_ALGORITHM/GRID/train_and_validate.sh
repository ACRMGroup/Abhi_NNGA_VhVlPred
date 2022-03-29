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


# ######################### SUB - ROUTINES SECTION ########################

do_set_path_environment_variables ()
{
   export LD_LIBRARY_PATH="/usr/local/lib:/acrm/gridengine/lib/grass6.0.0-i686-pc-linux-gnu-14_03_2005/lib/:/home/bsm2/abhi/lib:/home/bsm/martin/lib"

   export GENE_POPULATION=$genePopulation
   export TEMPORARY_SNNS_PATH="/acrm/home/abhi/tmp/SNNS/"
   export ABINTERFACE="/home/bsm2/abhi/VH-VL-INTERFACE/STRUCTURAL_ANALYSIS/"
   export ABHIDATADIR="/home/bsm2/abhi/DATA/"
   export PATH="/acrm/gridengine/bin/glinux:/usr/kerberos/bin:/acrm/usr/local/bin:/usr/bin:/usr/local/bin:/bin:/usr/X11R6/bin:/home/bsm2/abhi/BASIC_UTILITIES:/home/bsm/martin/bin:/home/bsm2/abhi/VH-VL-INTERFACE/KABAT_NUMBERING/"

}


do_set_other_environment_variables ()
{

   # Set number of bins if required

   if [ "$binOption" != "" ]									#  (1)
   then

      if [ "$numberOfBins" != "" ]								#  (2)
      then

         # Check if $numberOfBins is a valid number.

         expr $numberOfBins + 1 >& /dev/null

         if [ $? -ne 0 ]									#  (3)
         then

            # Here, the number of bins parameter is not a valid number. Hence, we export
            # a default value of 5.

            export NUMBER_OF_BINS=5

         else											# Else for (3)

            # The value in $numberOfBins is a valid number.

            export NUMBER_OF_BINS=$numberOfBins

         fi # End of "if [ $? -ne 0 ]"								# End of (3)

      else # End of "if [ "$numberOfBins != "" ]"						# Else for (2)

         # Here, $numberOfBins = ""
         # Set NUMBER_OF_BINS to $numberOfBins.

         export NUMBER_OF_BINS=$numberOfBins

      fi # End of "if [ "$numberOfBins" != "" ]"						# End of (2)

   fi # End of "if [ "$binOption" != "" ]"							# End of (1)

} # End of function "do_set_other_environment_variables".


# Parameters to the program:
#
# ./train_and_validate.sh $inputFilename $genePopulation $geneNumber $outputDirectory $generation $geneType

if [ $# -lt 7 ]
then
   echo
   echo "Usage: $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. Input filename"
   echo "2. Gene number (OR) Range (Eg. 1-200)"
   echo "3. Output directory"
   echo "4. Generation"
   echo "5. Gene type - P or C"
   echo "6. File with list of interface positions"
   echo "7. -rmse or -pear - Use Root Mean Square Error OR Pearsons correlation coefficient for calculation of score"
   echo "8. -binT or -binTV or -binVT or -binV - Use bins to represent the training and/or validation sets"
   echo "9. Number of bins (If using one of the bin options) - Default 5"
   echo
   exit 0
fi

inputFilename=$1
geneNumber=$2
outputDirectory=$3
generation=$4
geneType=$5
interfacePositionsListFilename=$6
scoreOption=$7
binOption=$8
numberOfBins=$9

if [ "$binOption" != "" ]
then
   if [ "$numberOfBins" == "" ]
   then
      echo
      echo "Please enter a value for number of bins"
      echo
      exit 0
   fi
fi

# Find number of genes in the input file.

genePopulation=`grep -c "GENE #" $inputFilename`

# Set the required environment variables.

do_set_path_environment_variables
do_set_other_environment_variables

# Check if the gene number parameter is proper.

expr $geneNumber + 1 >& /dev/null

if [ $? -ne 0 ]
then
   # In this case, $geneNumber is a range. Check if it is of the form 1-100.

   command=`echo $geneNumber | grep '-'`

   if [ "$command" == "" ]
   then
      echo
      echo "Invalid range: $geneNumber"
      exit 0
   fi

   # We now know that range is of the form 1-200. Check if the starting number
   # is less than the ending number. If not, report error and exit program.

   geneNumberStart=`echo $geneNumber | awk -F'-' '{print $1}'`
   geneNumberEnd=`echo $geneNumber | awk -F'-' '{print $2}'`

   if [ $geneNumberStart -gt $geneNumberEnd ]
   then
      echo
      echo "Invalid range $geneNumber as $geneNumberStart > $geneNumberEnd"
      echo
      exit 0
   fi

   # The gene string is a proper range.

   useBatch=1

else

   useBatch=0

fi

# Execute the program "train_and_validate.exe"

if [ "$useBatch" == "0" ]
then

   # In this case, score has to be calculated for a single gene.

   ./train_and_validate.exe -in $inputFilename -num $geneNumber -int $interfacePositionsListFilename -pdb $ABINTERFACE/Fv_Fab_pdb.lst -tor $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt -type $geneType -odir $outputDirectory -gen $generation $scoreOption $binOption

else

   # In this case, scores have to be calculated for a set of genes.

   geneNumber=$geneNumberStart

   if [ $geneNumberEnd -gt $genePopulation ]
   then
      geneNumberEnd=$genePopulation
   fi

   while [ $geneNumber -le $geneNumberEnd ]
   do

      ./train_and_validate.exe -in $inputFilename -num $geneNumber -int $interfacePositionsListFilename -pdb $ABINTERFACE/Fv_Fab_pdb.lst -tor $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt -type $geneType -odir $outputDirectory -gen $generation $scoreOption $binOption

      geneNumber=`expr $geneNumber + 1`

   done

fi

#! /bin/bash

# This program jumps starts the genetic algorithm with an initial set of
# parent genes and their scores.

# Parameters required for the program:
#

if [ $# -lt 4 ]
then
   echo
   echo "Usage: $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. Input file"
   echo "2. Clear existing text files - Y or N"
   echo "3. File number for creating an output file"
   echo "4. File with list of interface positions"
   echo "5. -pear OR -rmse - procedure for calculating score of a gene"
   echo "6. -binT or -binTV or -binV or -binVT - Optional"
   echo "7. Number of bins (If the bin option is used) - Optional"
   echo
   exit 0
fi

inputFilename=$1
clearTextFilesOption=$2
fileNumber=$3
interfacePositionsListFilename=$4
scoringOption=$5
binOption=$6
numberOfBins=$7

echo "NUmber of bins: $numberOfBins"
read

# Check if the number of bins 

if [ "$binOption" != "" ]
then
   if [ "$numberOfBins" == "" ]
   then
      echo
      echo "Please input number of bins"
      echo
      exit 0
   fi
fi 

case "$binOption" in
   "-binV") ;;
   "-binTV") ;;
   "-binV") ;;
   "-binVT") ;;
   "-binT") ;;
   "") numberOfBins=0 ;;
   *) echo
      echo "Invalid option \"$binOption\""
      echo
      exit 0
      ;;
esac

userName=`whoami`

if [ "$clearTextFilesOption" == "Y" ]
then
   clearTextFiles="Y"
else
   clearTextFiles="N"
fi

# Set environment variables and the output directory.

genePopulation=`grep -c "GENE #" $inputFilename`

export GENE_POPULATION=$genePopulation

outputDirectory="./OUTPUTS"

# Print all the parameters.

echo
echo "Program: $0"
echo "Input file: $inputFilename"
echo "Clear text files option: $clearTextFilesOption"
echo "Number of genes in file \"$inputFilename\": $genePopulation"
echo "-----------------------"

# Check if input file is present

if [ ! -e $inputFilename ]
then
   echo
   echo "File \"$inputFilename\" does not exist. Aborting program"
   echo
   exit 0
fi

# Check if the environment variable TEMPORARY_SNNS_PATH is set. If it isn't, then
# report an error message and exit from program.

if [ ! `echo $TEMPORARY_SNNS_PATH` ]
then
   echo
   echo "Environment variable \"TEMPORARY_SNNS_PATH\" not set. Aborting program"
   echo
   exit 0
fi

# Print the parameters.

# Clear existing files.

if [ "$clearTextFiles" == "Y" ]
then
   mv $jumpStartFilename $$.tmp
   rm -f *txt
   mv $$.tmp $jumpStartFilename
fi

find $outputDirectory | xargs rm -f

# Set alias for cp command

alias cp=cp

# Remove unnecessary files created by gridengine

find ./ -name "train_and_validate.sh.*" | xargs rm -f | sed 's/ //g'

find $TEMPORARY_SNNS_PATH | xargs rm -f

# Evaluate the genes in the input file.

generation=1
geneType=P

genePopulation=`grep -c "GENE #" $inputFilename`

geneNumber=1

# ./train_and_validate.sh $inputFilename $genePopulation $geneNumber $outputDirectory $generation $geneType

while [ $geneNumber -le $genePopulation ]
do
   if [ "$numberOfBins" != "0" ]
   then

   qsub -cwd -l s_rt=00:02:00 ./train_and_validate.sh $inputFilename $geneNumber $outputDirectory $generation $geneType $interfacePositionsListFilename $scoringOption $binOption $numberOfBins

   else

   qsub -cwd -l s_rt=00:02:00 ./train_and_validate.sh $inputFilename $geneNumber $outputDirectory $generation $geneType $interfacePositionsListFilename $scoringOption

   fi

   geneNumber=`expr $geneNumber + 1`
done

# Monitor the grid queue till calculation of scores of all child genes is complete.

temp="dummy"
while [ 1 ]
   [ "$temp" != "" ]
do
   temp=`qstat | grep $userName`
done

# Check if number of output files equals the gene population

numberOfOutputFiles=`find $outputDirectory -name "P*.txt" | wc -l`

if [ $numberOfOutputFiles -ne $genePopulation ]
then
   echo
   echo
   echo "Number of output files does not match number of genes"
   echo
   echo
   exit 0
fi

# Remove files that arent needed.

find ./ -name "train_and_validate.sh.*" | xargs rm -f

find $TEMPORARY_SNNS_PATH | xargs rm -f

# Concatenate all the genes and their scores into one file.

if [ "$binOption" != "" ]
then
   outputFilename="genes_and_scores_"$fileNumber$binOption".txt"
else
   outputFilename="genes_and_scores_"$fileNumber".txt"
fi

rm -f $outputFilename

for file in `find $outputDirectory -name "P*.txt"`
do
   cat $file >> $outputFilename
   rm -f $file
done

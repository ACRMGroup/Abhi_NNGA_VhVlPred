#! /bin/bash

# This program jumps starts the genetic algorithm with an initial set of
# parent genes and their scores.

# Parameters required for the program:
#

if [ $# -lt 5 ]
then
   echo
   echo "Usage: $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. Input file"
   echo "2. Selection procedure - Roulette or Rank or int"
   echo "3. Mutation rate"
   echo "4. Clear existing text files - Y or N"
   echo "5. File number for creating an output file"
   echo "6. Gene number to start with"
   echo "7. Number of bins - Optional"
   echo
   exit 0
fi

inputFilename=$1
selectionProcedure=$2
mutationRate=$3
clearTextFilesOption=$4
fileNumber=$5
geneStartNumber=$6

if [ -n $7 ]
then
   numberOfBins=$7
fi

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
export MUTATION_RATE=$mutationRate

outputDirectory="./OUTPUTS"

# Print all the parameters.

echo
echo "Program: $0"
echo "Input file: $inputFilename"
echo "Selection procedure; $selectionProcedure"
echo "Mutation rate: $mutationRate"
echo "Clear text files option: $clearTextFilesOption"
echo "Number of genes in file \"$inputFilename\": $genePopulation"
echo "Gene start number: $geneStartNumber"
echo "Number of bins: $numberOfBins"
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

if [ -n $geneStartNumber ]
then
   geneNumber=$geneStartNumber
else
   geneNumber=1
fi

# ./train_and_validate.sh $inputFilename $genePopulation $geneNumber $outputDirectory $generation $geneType $mutationRate

while [ $geneNumber -le $genePopulation ]
do
   if [ -n $numberOfBins ]
   then

      qsub -cwd -p -1023 train_and_validate.sh $inputFilename $genePopulation $geneNumber $outputDirectory $generation $geneType $mutationRate $numberOfBins

   else

      qsub -cwd -p -1023 train_and_validate.sh $inputFilename $genePopulation $geneNumber $outputDirectory $generation $geneType $mutationRate

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

outputFilename="genes_and_scores_"$fileNumber".txt"

for file in `find $outputDirectory -name "P*.txt"`
do
   cat $file >> $outputFilename
   rm -f $file
done

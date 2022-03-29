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
   echo "1. Number of genes"
   echo "2. Number of generations"
   echo "3. Selection procedure - Roulette or Rank"
   echo "4. File to jump start with"
   echo "5. Mutation rate"
   echo "6. Clear existing text files - Y or N (Optional parameter)"
   echo
   exit 0
fi

genePopulation=$1
numberOfGenerations=$2
selectionProcedure=$3
jumpStartFilename=$4
mutationRate=$5

if [ "$6" == "Y" ]
then
   clearTextFiles="Y"
else
   clearTextFiles="N"
fi

export GENE_POPULATION=$genePopulation
export MUTATION_RATE=$mutationRate

outputDirectory="./OUTPUTS"

# Check if selectionProcedure matches either "Roulette" or "Rank"

# if [ "$selectionProcedure" != "Rank" ]
# then
#    if [ "$selectionProcedure" != "Roulette" ]
#    then
#       echo
#       echo "Invalid parameter \"$selectionProcedure\" for selection procedure"
#       echo
#       exit 0
#    fi
# fi

# Check if jump start file is present

temp=`ls -1 $jumpStartFilename`

if [ "$temp" == "" ]
then
   echo
   echo "File \"$jumpStartFilename\" does not exist. Aborting program"
   echo
   exit 0
fi

# Print the parameters.

echo "Program: $0"
echo "Number of genes: $genePopulation"
echo "Number of generations: $numberOfGenerations"
echo "Selection procedure: $selectionProcedure"
echo "Jump start file: $jumpStartFilename"
echo "Mutation rate: $mutationRate"
echo "Clear text files: $clearTextFiles"
echo "----------------------"
echo
echo

# Clear existing files.

if [ "$clearTextFiles" == "Y" ]
then
   mv $jumpStartFilename $$.tmp
   rm -f *txt
   mv $$.tmp $jumpStartFilename
fi

rm -f $outputDirectory/*

# Set alias for cp command

alias cp=cp

# Remove unnecessary files created by gridengine

rm -f train_and_validate.sh.e*
rm -f train_and_validate.sh.o*

rm -rf $TEMPORARY_SNNS_PATH
mkdir $TEMPORARY_SNNS_PATH

# Create child genes over N generations.

generation=1
geneType=C

while [ $generation -le $numberOfGenerations ]
do
   # First, create the child genes.

   ./create_child_genes.exe -in $jumpStartFilename -gen $generation -proc $selectionProcedure

   # Calculate the scores of the child genes.

   geneNumber=1

   while [ $geneNumber -le $genePopulation ]
   do
      qsub -cwd train_and_validate.sh child_genes.txt $genePopulation $geneNumber $outputDirectory $generation $geneType
      geneNumber=`expr $geneNumber + 1`
   done

   # Monitor the grid queue till calculation of scores of all child genes is complete.

   temp="dummy"
   while [ 1 ]
      [ "$temp" != "" ]
   do
      temp=`qstat | grep abhi`
   done

   # Remove files that arent needed.

   rm -f train_and_validate.sh.e*
   rm -f train_and_validate.sh.o*

   rm -rf $TEMPORARY_SNNS_PATH
   mkdir $TEMPORARY_SNNS_PATH/

   # Concatenate all the child genes into a single file.

   cat $outputDirectory/C*.txt > child_genes_and_scores.txt

   # Write the child genes into another file.

   echo "----------------" >> child_genes_and_scores_copy.txt
   echo "GENERATION #$generation" >> child_genes_and_scores_copy.txt
   echo "----------------" >> child_genes_and_scores_copy.txt
   echo >> child_genes_and_scores_copy.txt

   cat $outputDirectory/C*.txt >> child_genes_and_scores_copy.txt

   echo >> child_genes_and_scores_copy.txt

   # Select the best genes from the pool of parent and child genes.

   ./select_best_genes.exe -mom $jumpStartFilename -kid child_genes_and_scores.txt -out $jumpStartFilename -gen $generation

   # Copy parent genes and their scores into the file "parent_genes_and_scores_copy.txt"

   echo "----------------" >> parent_genes_and_scores_copy.txt
   echo "GENERATION #$generation" >> parent_genes_and_scores_copy.txt
   echo "----------------" >> parent_genes_and_scores_copy.txt
   echo >> parent_genes_and_scores_copy.txt
   cat parent_genes_and_scores.txt >> parent_genes_and_scores_copy.txt

   # Increment generation and iterate the process till N generations are complete.

   generation=`expr $generation + 1`

done

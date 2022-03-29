#! /bin/bash

# This script wraps around programs that implement a genetic algorithm. The steps involved in the program
# are as follows:

# *) Create parent genes using the program create_parent_genes.exe. Genes are written to file "parent_genes.txt".

# *) Calculate scores for the parent genes using the program "train_and_validate.exe". This program
#    is run for every gene in the file "parent_genes.txt". The execution of the program on each gene
#    is carried out in parallel on the grid.
#
# *) Create child genes from the file "parent_genes_and_scores.txt" using the program "create_child_genes.exe".
#    Child genes are written to the file "child_genes.txt". The score of each gene is written to a seperate file.
#
# *) Calculate scores for the child genes using the program "train_and_validate.exe". This program
#    is run for every gene in the file "child_genes.txt". The execution of the program on each gene
#    is carried out in parallel on the grid.
#
# *) Compare the scores of parent genes and child genes and select a set of genes that will become the
#    parent genes for the next generation. The selection of genes can be based on one of several approaches.
#    The old genes in the file "parent_genes_and_scores.txt" are replaced by the new set of parent genes.
#
# *) Iterate steps 2 to 4 for N generations of the execution.

# Parameters required for the program:
#
# Number of genes to be created.
# Number of generations over which the GA has to be run
# Selection procedure while creating the child genes - Roulette wheel or Rank based.
# Mutation rate - Probability of making a mutation in every allele in a gene.

if [ $# -lt 4 ]
then
   echo
   echo "Usage: $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. Number of genes"
   echo "2. Number of generations"
   echo "3. Selection procedure - Roulette or Rank or int"
   echo "4. Mutation rate"
   echo
   exit 0
fi

genePopulation=$1
numberOfGenerations=$2
selectionProcedure=$3
mutationRate=$4

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
#    else
#    if [ "$selectionProcedure" != "int" ]
#    then
#       echo
#       echo "Invalid parameter \"$selectionProcedure\" for selection procedure"
#       echo
#       exit 0
#    fi
# fi

# Print the parameters.

echo "Program: $0"
echo "Number of genes: $genePopulation"
echo "Number of generations: $numberOfGenerations"
echo "Selection procedure: $selectionProcedure"
echo "Mutation rate: $mutationRate"
echo "----------------------"
echo
echo

# Clear existing files.

rm -f *.txt $outputDirectory/*

# Set alias for cp command

alias cp=cp

# Create the parent genes.

./create_parent_genes.exe -int ../all_interface_positions.dat

# Calculate scores for parent genes.

geneNumber=1
geneType=P
generation=1

while [ $geneNumber -le $genePopulation ]; do

   qsub -cwd train_and_validate.sh parent_genes.txt $genePopulation $geneNumber $outputDirectory $generation $geneType

   geneNumber=`expr $geneNumber + 1`

done

# Monitor the grid queue using qstat to ensure that all processes have been completed before
# proceeding to the next step.

temp="dummy"
while [ 1 ]
   [ "$temp" != "" ]
do
   temp=`qstat | grep abhi`
done

# Check number of output files created in Output directory.

numberOfOutputFiles=`ls -1 $outputDirectory/P* | wc -l`

if [ "$numberOfOutputFiles" != "$genePopulation" ]
then
   echo
   echo
   echo "Number of output files does not match number of genes"
   echo
   echo
   exit 0
fi

# Remove unnecessary files created by gridengine

rm -f train_and_validate.sh.e*
rm -f train_and_validate.sh.o*

for file in `ls -1 $TEMPORARY_SNNS_PATH/`
do
   rm -f $file
done

# Concatenate all the parent scores into a single file - parent_genes_and_scores.txt

cat $outputDirectory/P*.txt > parent_genes_and_scores.txt

# Create child genes over N generations.

generation=1
geneType=C

while [ $generation -le $numberOfGenerations ]
do
   # First, create the child genes.

   ./create_child_genes.exe -in parent_genes_and_scores.txt -gen $generation -proc $selectionProcedure

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
   mkdir $TEMPORARY_SNNS_PATH

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

   ./select_best_genes.exe -mom parent_genes_and_scores.txt -kid child_genes_and_scores.txt -out parent_genes_and_scores.txt -gen $generation

   # Copy contents of parent_genes_and_scores.txt into file parent_genes_and_scores_copy.txt

   echo "----------------" >> parent_genes_and_scores_copy.txt
   echo "GENERATION #$generation" >> parent_genes_and_scores_copy.txt
   echo "----------------" >> parent_genes_and_scores_copy.txt
   echo >> parent_genes_and_scores_copy.txt
   cat parent_genes_and_scores.txt >> parent_genes_and_scores_copy.txt

   # Increment generation and iterate the process till N generations are complete.
 
   generation=`expr $generation + 1`

done

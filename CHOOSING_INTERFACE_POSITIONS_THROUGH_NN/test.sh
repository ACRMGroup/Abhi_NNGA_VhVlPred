#! /bin/sh

if [ $# -lt 3 ]
then
   echo
   echo "Usage: $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. Name of file with the interface positions"
   echo "2. String for training set boundary (Eg. 126-124)"
   echo "3. File with numbering for all antibody structures"
   echo "4. Output file (Optional) - Default: errors.txt in current directory"
   echo "5. Name of untrained network file (Optional) - Default: with 10 hidden nodes"
   echo "6. Number of cycles for training (Optional) - Default: 200"
   echo
   exit 0
fi

interfacePositionsListFilename=$1
trainingSetBoundaryString=$2
allAntibodiesNumberingFilename=$3

if [ $# -gt 3 ]
then
   outputFilename=$4
else
   outputFilename="errors.txt"
fi

if [ $# -gt 4 ]
then
   untrainedNetworkFilename=$5
else
   untrainedNetworkFilename="UNTRAINED_NEURAL_NETWORKS/1i_10h.net"
fi

if [ $# -gt 5 ]
then
   numberOfCycles=$6
else
   numberOfCycles=200
fi

if [ ! -e $interfacePositionsListFilename ]
then
   echo
   echo "File $interfacePositionsListFilename does not exist"
   echo
   exit 0
fi

# Come to the actual part of doing the job.

temporaryInterfaceFilename="TEMPORARY_FILES/interface.lst"

rm -f $outputFilename

echo "Number of cycles: $numberOfCycles" >> $outputFilename
echo "Untrined network file: $untrainedNetworkFilename" >> $outputFilename
echo "------------------------------------" >> $outputFilename
echo >> $outputFilename
echo >> $outputFilename

for interfacePosition in `cat $interfacePositionsListFilename`
do
   echo $interfacePosition

   c=`grep -c "$interfacePosition " $allAntibodiesNumberingFilename`

   echo "Interface position: $interfacePosition" >> $outputFilename

   if [ "$c" = "481" ]
   then

      echo "-----------------------------"
      echo "Processing $interfacePosition"
      echo "-----------------------------"

      echo $interfacePosition > $temporaryInterfaceFilename

      intermediateFilename="INTERMEDIATE_FILES/$interfacePosition"_intermediate.txt

      sh create_intermediate_file.sh TEMPORARY_FILES/interface.lst $intermediateFilename

      cmd="./final.exe -int $intermediateFilename -net $untrainedNetworkFilename -trb $trainingSetBoundaryString -out $outputFilename -cyc $numberOfCycles"

      `$cmd`

      echo "##########################" >> $outputFilename

      rm -f $intermediateFilename

   else
      echo "Only $c numbering entries found for position $interfacePosition" >> $outputFilename
      echo "Unable to process further" >> $outputFilename
      echo "##########################" >> $outputFilename
   fi

done

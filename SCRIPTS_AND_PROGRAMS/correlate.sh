if [ $# -lt 1 ]
then
   echo
   echo "Usage: sh $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "a) 10 - Network with 10 hidden nodes"
   echo "b) 20 - Network with 20 hidden nodes"
   echo "c) 30 - Network with 30 hidden nodes"
   echo "4) 40 - All the above"
   echo
   exit 0
fi

baseDirectory=$HOME/SNNS

input=$1


for validationPatternFilename in `ls -1 $HOME/SNNS/PATTERN_FILES/*_validation.pat`
do
   # If pattern file is frequency_average_validation.pat, then we need the prefix "frequency_average".

   prefix1=`basename $validationPatternFilename .pat`	# prefix - frequency_average
   prefix=`echo $prefix1 | sed 's/_validation//'`

   # 1. 20 Input nodes and 10 Hidden nodes.

   if [ $input -gt 9 ]
   then

      networkType="20i_10h"
   
      neuralNetworkFilename="$baseDirectory/TRAINED_NEURAL_NETWORKS/$prefix"_"trained"_"$networkType".net
      neuralNetworkPredictionsOutputFilename="$baseDirectory/NEURAL_NETWORK_PREDICTIONS/$prefix"_"$networkType".res
      angleCorrelationFilename="$baseDirectory/ANGLE_COMPARISON_RESULT_FILES/$prefix"_"$networkType".out

      perl $baseDirectory/SCRIPTS_AND_PROGRAMS/validate.pl $neuralNetworkFilename $validationPatternFilename $neuralNetworkPredictionsOutputFilename $angleCorrelationFilename

   fi

   # 2. 20 Input nodes and 20 Hidden nodes.

   if [ $input -gt 19 ]
   then

      networkType="20i_20h"

      neuralNetworkFilename="$baseDirectory/TRAINED_NEURAL_NETWORKS/$prefix"_"trained"_"$networkType".net
      neuralNetworkPredictionsOutputFilename="$baseDirectory/NEURAL_NETWORK_PREDICTIONS/$prefix"_"$networkType".res
      angleCorrelationFilename="$baseDirectory/ANGLE_COMPARISON_RESULT_FILES/$prefix"_"$networkType".out

      perl $baseDirectory/SCRIPTS_AND_PROGRAMS/validate.pl $neuralNetworkFilename $validationPatternFilename $neuralNetworkPredictionsOutputFilename $angleCorrelationFilename

   fi

   # 3. 20 Input nodes and 30 Hidden nodes.

   if [ $input -gt 29 ]
   then

      networkType="20i_30h"

      neuralNetworkFilename="$baseDirectory/TRAINED_NEURAL_NETWORKS/$prefix"_"trained"_"$networkType".net
      neuralNetworkPredictionsOutputFilename="$baseDirectory/NEURAL_NETWORK_PREDICTIONS/$prefix"_"$networkType".res
      angleCorrelationFilename="$baseDirectory/ANGLE_COMPARISON_RESULT_FILES/$prefix"_"$networkType".out

      perl $baseDirectory/SCRIPTS_AND_PROGRAMS/validate.pl $neuralNetworkFilename $validationPatternFilename $neuralNetworkPredictionsOutputFilename $angleCorrelationFilename

   fi

done

if [ $# -lt 1 ]
then
   echo "Usage: sh $0 <Name of C program> <Name of executable (Optional)>"
   exit 0
fi

inputCfile=$1
prefix=`basename $inputCfile .c`

if [ "$2" != "" ]
then
   executable=$2
else
   executable=$prefix.exe
fi

gcc -ansi -Wall -pedantic -g $inputCfile  -o $executable -lgsl -lgslcblas -L ~/lib/ -labs -I ~/include/ -L ~martin/lib/ -lbiop -lgen -I ~martin/include -L /LINUX/local/grass6.0.0-i686-pc-linux-gnu-14_03_2005/lib/ -lgrass -I /LINUX/local/grass6.0.0-i686-pc-linux-gnu-14_03_2005/include/

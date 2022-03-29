if [ $# -lt 1 ]
then
   echo "Please enter name of the input C file"
   exit 0
fi

inputCfile=$1
prefix=`basename $inputCfile .c`
executable=$prefix.exe

gcc -ansi -Wall -pedantic -g $inputCfile  -o $executable -L ~/lib/ -labs -I ~/include/ -L ~martin/lib/ -lbiop -lgen -I ~martin/include -L /LINUX/local/grass6.0.0-i686-pc-linux-gnu-14_03_2005/lib/ -lgrass -I /LINUX/local/grass6.0.0-i686-pc-linux-gnu-14_03_2005/include/

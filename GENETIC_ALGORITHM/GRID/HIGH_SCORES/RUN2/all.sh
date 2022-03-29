#! /bin/sh

sh bruteforce.sh best_1.txt int 0.001 N 1 2 >& 1.log

find OUTPUTS/ | xargs rm -f

sh bruteforce.sh best_2.txt int 0.001 N 2 1 >& 2.log

find OUTPUTS/ | xargs rm -f

sh bruteforce.sh best_3.txt int 0.001 N 3 1 >& 3.log

find OUTPUTS/ | xargs rm -f

sh bruteforce.sh best_4.txt int 0.001 N 4 1 >& 4.log

find OUTPUTS/ | xargs rm -f

sh bruteforce.sh best_5.txt int 0.001 N 5 1 >& 5.log

find OUTPUTS/ | xargs rm -f

sh bruteforce.sh best_6.txt int 0.001 N 6 1 >& 6.log

find OUTPUTS/ | xargs rm -f

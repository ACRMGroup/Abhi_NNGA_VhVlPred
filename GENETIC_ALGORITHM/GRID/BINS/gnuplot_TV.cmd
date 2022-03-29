set terminal png
set out 'TV.png'

set title "Bin style representation: TRAINING AND VALIDATION PATTERNS"
set xlabel "Correlation coefficient after using bins"
set ylabel "Correlation coefficient without using bins"

plot [0.5:0.7] [0.64:0.69] 'TV_5_bin_actual.txt' using 1:2 w points title "5 bins",'TV_10_bin_actual.txt' using 1:2 w points title "10 bins",'TV_15_bin_actual.txt' using 1:2 w points title "15 bins",'TV_20_bin_actual.txt' using 1:2 w points title "20 bins",'TV_25_bin_actual.txt' using 1:2 w points title "25 bins"

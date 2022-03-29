set terminal png
set out 'V.png'

set title "Bin style representation: VALIDATION PATTERNS"
set xlabel "Correlation coefficient after using bins"
set ylabel "Correlation coefficient without using bins"

plot [0.5:0.7] [0.64:0.69] 'V_5_bin_actual.txt' using 1:2 w points title "5 bins",'V_10_bin_actual.txt' using 1:2 w points title "10 bins",'V_15_bin_actual.txt' using 1:2 w points title "15 bins",'V_20_bin_actual.txt' using 1:2 w points title "20 bins",'V_25_bin_actual.txt' using 1:2 w points title "25 bins"

set terminal pngcairo size 1280,720 enhanced font 'Arial,12'
set output "residuals.png"

set datafile separator ","
set title "Konvergence výpočtu – reziduum"
set xlabel "Iterace"
set ylabel "Reziduum (L∞)"
set grid
set logscale y

plot "residuals.csv" using 1:2 with lines lw 2 title "reziduum"

set output
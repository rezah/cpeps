set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "1/N"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"
#set ytics nomirror (10.3, 10.4,10.5,10.6, 10.7, 10.8, 10.9);


set xtics ( "1/4" 0.25,"1/8" 0.12, "1/12" 0.08);
set tics scale 3

set key box t c 
set key font ",30"
set key spacing 1.30
set key width 0.90


set output "Density.eps"



p [0.05:0.29]  [-15:-14.5] -14.90 t "QMC" with line lw 4 linetype 7 lc rgb "#0251fa", "N4.txt" u ($1):($5)  t  "N_e=4" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#ed09e2" ,"N10.txt" u ($1):($5)  t  "N_e=10" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ed4a09","N20.txt" u ($1):($5)  t  "N_e=20" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#09a5ed" 



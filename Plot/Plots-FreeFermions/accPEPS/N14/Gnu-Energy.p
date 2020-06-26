set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "iteration"   font ",30" textcolor rgb "black"
set ylabel ""  font ",30"  textcolor rgb "black"
set ytics nomirror (45,46,47,48,57,60);


set xtics ( 5,10,20,30,40);
set tics scale 3

set key box b l at 9,48.
set key font ",30"
set key spacing 1.30
set key width 0.90


set output "Density.eps"



p [0.5:41]  [43.9:49.] 44.936381 notitle with line lw 4 linetype 7 lc rgb "#0251fa", "FU.txt" u ($1):($2)  t  "{/Symbol c}_o=8, D=8" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#ed09e2" , "FUPEPS.txt" u ($1):($2)  t  "{/Symbol c}_b=150, D=6" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#d91111", "FUPEPSD8.txt" u ($1):($2)  t  "{/Symbol c}_b=250, D=8" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#0999db"


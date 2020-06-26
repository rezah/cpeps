set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "iteration"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"
#set ytics nomirror (10.3, 10.4,10.5,10.6, 10.7, 10.8, 10.9);


set xtics ( 40, 80, 160,260);
set tics scale 3

set key box t c
set key font ",30"
set key spacing 1.30
set key width 0.90


set output "Density.eps"



p [1:340]  "EnergyRGv1.txt" u ($1):($3)  t  "D=4" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#ed09e2" ,"EnergyRGv2.txt" u ($1+220):($3+0.004)  t "D=6" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#ed4a09"

#,"N20eta1.txt" u ($1):($3)  t  "{/Symbol h}=1,N_e=20" with points  pointtype 19 lw 3  ps 4.0 lc rgb "#ed09e2" 



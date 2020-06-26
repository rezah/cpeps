set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "{/Symbol h}"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"
#set ytics nomirror (10.3, 10.4,10.5,10.6, 10.7, 10.8, 10.9);


set xtics ( -1.5, 0, 1.5);
set tics scale 3

set key box c c
set key font ",30"
set key spacing 1.30
set key width 0.90


set output "Density.eps"



p [-1.5:4]  "Eqmc.txt" u ($1):($2):($3)  t  "QMC" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#ed09e2" ,"Epeps.txt" u ($1):($2)  t  "TN" with points  pointtype 8 lw 3  ps 4.0 lc rgb "#09a1ed"


#,"N20eta1.txt" u ($1):($3)  t  "{/Symbol h}=1,N_e=20" with points  pointtype 19 lw 3  ps 4.0 lc rgb "#ed09e2" 



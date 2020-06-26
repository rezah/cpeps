set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "1/{/Symbol r}"   font ",30" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",30"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics ( 3,2, 1,0.5);
set tics scale 4

set key box at 1.70,-4.10
set key font ",30"
set key spacing 1.05
set key width 1.


set output "Density.eps"


set multiplot

p [0.3:3.6] [ -4.5 : -0.80] "FU22.txt" u (6./2.):((log(($2-$3)*(1/$3))/log(10)))  t  "{/Symbol e}=0.66" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4" ,"FU33.txt" u (6./3.):((log(($2-$3)*(1/$3))/log(10)))  notitle with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4","FU44.txt" u (6./4.):((log(($2-$3)*(1/$3))/log(10)))  notitle with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4","FU66.txt" u (6./6.):((log(($2-$3)*(1/$3))/log(10)))  notitle with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4","FU88.txt" u (6./8.):((log(($2-$3)*(1/$3))/log(10)))  notitle with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4","FU100.txt" u (6./10.):((log(($2-$3)*(1/$3))/log(10)))  notitle with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4","FU144.txt" u (6./14.):((log(($2-$3)*(1/$3))/log(10)))  notitle with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4"


set origin .39, .3750
set size 0.55,0.550
clear
unset key
unset xlabel
set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "1/{/Symbol c}_{o}"   font ",30" textcolor rgb "black"
set ylabel ""  font ",30"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);

set xtics ( "1/4" 0.25,"1/3" 0.333, "1/8" 0.125, "1/16" 0.0625)
set tics scale 3

set key box at  0.3,-2.5
set key font ",25"
set key spacing 1.05
set key width 1.02


p [0.02:0.31] [ -3.50 : -0.1] "FU4.txt" u (1/$1):((log(($2-9.950227785)*(1/9.950227785))/log(10)))  t "{/Symbol e}=0.86" with points  pointtype 13 lw 3  ps 2.50 lc rgb "#204a87","FU2.txt" u (1/$1):((log(($2-10.895071280)*(1/10.895071280))/log(10)))  t "{/Symbol e}=0.66" with points  pointtype 11 lw 3  ps 2.50 lc rgb "#f57900","FU6.txt" u (1/$1):((log(($2-10.9491066261)*(1/10.9491066261))/log(10)))  t "{/Symbol e}=0.46" with points  pointtype 15 lw 3  ps 2.50 lc rgb "#f53500" 

unset multiplot


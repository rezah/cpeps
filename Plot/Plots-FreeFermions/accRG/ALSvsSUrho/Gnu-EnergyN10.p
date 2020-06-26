set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/{/Symbol c}_{r}"   font ",40" textcolor rgb "black"
set ylabel ""  font ",40"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics ("1/2" 0.5,"1/4" 0.25, "1/6" 0.166,"1/10" 0.1);
set tics scale 4

set key box b r 
set key font ",35"
set key spacing 1.2
set key width 1.20


set output "RGN10.eps"



p [ 0.05: 0.54] [ -4.0 : -0.1] "accenv.txt" u ($1):((log(($3+10.8068766146)*(1/10.8068766146))/log(10)))  t "SVD" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#3465a4" ,"accenv.txt" u ($1):((log(($2+10.8068766146)*(1/10.8068766146))/log(10)))  t "ALS" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#a4004d" ,"rand10.txt" u ($1):((log(($2+607.74274)*(1/607.74274))/log(10)))  t "rand" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ed18a7" 




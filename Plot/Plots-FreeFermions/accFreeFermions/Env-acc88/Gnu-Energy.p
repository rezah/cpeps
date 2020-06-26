set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/{/Symbol c}"   font ",35" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",35"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics ( "1/4" 0.25,"1/6" 0.1666, "1/8" 0.125, "1/5" 0.2)
set tics scale 4

set key box b r 
set key font ",30"
set key spacing 1.2
set key width 1.20


set output "Density.eps"



p [0.1:0.3] [ -3.80 : -1.50] "FU4.txt" u (1/$1):((log(($2+0.75356677360)*(1/0.7535667))/log(10)))  t "4{/Symbol \264}4" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#f57900"


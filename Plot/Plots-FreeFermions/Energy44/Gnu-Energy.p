set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "iteration"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",40"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);



set xtics (100,200,300);
set tics scale 4

set key box b l 
set key font ",30"
set key spacing 1.5
#set key width 2.5


set output "Energy.eps"



p [1 : 250] [ -3.3 : -0.5] -1.3 t "{/Symbol c}_{b}=12" with line lw 4 linetype 7 lc rgb "#006400","FUD5.txt" u ($1):((log(($3-8.973011)*(1/8.973011))/log(10)))  t "D=5" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#5c3566" ,"FUD7.txt" u (($1)+100):((log(($3-8.973011)*(1/8.973011))/log(10)))  t "D=7" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#21908d","FUD8.txt" u (($1)+210):((log(($3-9.005173011)*(1/9.005173011))/log(10)))  t "D=12" with points  pointtype 5 lw 3  ps 4.0 lc rgb "#f57900"



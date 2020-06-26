set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "iteration"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",40"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);



set xtics (0,20,40,60,80,100);
set tics scale 4

set key box b l 
set key font ",30"
set key spacing 1.5
#set key width 2.2


set output "Energy.eps"



p [1 : 60] [ -4 : -0.5] "FUD4.txt" u ($1):((log(($2-897.30)*(1/897.30))/log(10)))  t "D=4" with points  pointtype 5 lw 3  ps 4.0 lc rgb "#f57900" ,"FUD5.txt" u (($1)+10):((log(($2-897.30)*(1/897.30))/log(10)))  t "D=5" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#5c3566","FUD6.txt" u (($1)+21):((log(($2-897.30)*(1/897.30))/log(10)))  t "D=6" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#cc0000" ,"FUD7.txt" u (($1)+39):((log(($3-897.30)*(1/897.30))/log(10)))  t "D=7" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#3465a4",-1.18 t "p=8, D=4" with line lw 4 linetype 7 lc rgb "#006400",-1.2 t "p=12, D=6" with line lw 4 linetype 7 lc rgb "#3465a4"



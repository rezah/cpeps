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



p [1 : 250] [ -4 : -0.5] -1.3 t "{/Symbol c}_{b}=12" with line lw 4 linetype 7 lc rgb "#006400", "FUD4.txt" u ($1):((log(($2-7.263932)*(1/7.263932))/log(10)))  t "D=4" with points  pointtype 5 lw 3  ps 4.0 lc rgb "#f57900" ,"FUD5.txt" u (($1)+50):((log(($2-7.263932)*(1/7.263932))/log(10)))  t "D=5" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#5c3566","FUD6.txt" u (($1)+128):((log(($2-7.263932)*(1/7.263932))/log(10)))  t "D=6" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#cc0000" ,"FUD7.txt" u (($1)+200):((log(($2-7.263932)*(1/7.263932))/log(10)))  t "D=7" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#21908d"



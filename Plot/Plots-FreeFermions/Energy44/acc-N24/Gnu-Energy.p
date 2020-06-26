set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "D"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",40"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);



#set xtics (0,20,40,60,80,100);
set tics scale 4

set key box b l 
set key font ",30"
set key spacing 1.5
#set key width 2.2


set output "Energy.eps"



p [3 : 16] [ -3.8 : -0.5] "U1.txt" u ($1):((log(($2-5.0278640)*(1/5.0278640))/log(10)))  t "N=4*4" with points  pointtype 5 lw 3  ps 4.0 lc rgb "#f57900" ,"Z2.txt" u ($1):((log(($2-5.0278640)*(1/5.0278640))/log(10)))  t "N=8*8" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#5c3566" 




set terminal postscript eps size 6.8,6.8 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "x"   font ",35" textcolor rgb "black"
set ylabel "y"  font ",35"  textcolor rgb "black"
set zlabel "<n_i>"  font ",35"  textcolor rgb "black"
set xtics (0.2,0.4,0.6,0.8,2.0,2.6);
set ytics (0.2,0.4,0.6,0.8,2.0,2.6);
#set ztics (0.0, 0.05,0.1,0.15,0.2,0.25,0.3);
set ztics ( 0.02 ,0.04,0.06,0.08,0.1,0.12);



set tics scale 2

set key box t l 
set key font ",30"
set key spacing 1.5
set key width 2.2


#set output "particle4.eps"
set output "particle8.eps"
#set output "particle16.eps"


splot "particle8.txt" u 1:2:3   t "PEPS" with points palette pointsize 3 pointtype 7, 0.05*sin(x*pi)*sin(y*pi) t "sin(x)*sin(y)" with lines

#splot "particle4.txt" u 1:2:3   t "PEPS" with points palette pointsize 3 pointtype 7, 0.15*sin(x*pi)*sin(y*pi) t "sin(x)*sin(y)" with lines

#{/Symbol @\326\140}2


set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "x"   font ",35" textcolor rgb "black"
set ylabel "y"  font ",35"  textcolor rgb "black"
set zlabel "<n_i>"  font ",35"  textcolor rgb "black"
set xtics (0.4,0.8,1.2,1.6,2,2.4,2.6,3);
set ytics (0.4,0.8,1.2,1.6,2,2.4,2.6,3);
#set ztics (0.0, 0.05,0.1,0.15,0.2,0.25,0.3);
set ztics ( 0.02 ,0.04,0.06,0.08,0.1,0.12);



set tics scale 2

set key box t l 
set key font ",30"
set key spacing 1.5
set key width 2.2


#set output "particle10.eps"
#set output "particle8.eps"
set output "Potential10.eps"

set xrange [0.2:4]
set yrange [0.2:4]
#set zrange [0:0.2]

set view map
splot "Potential.txt" u 1:2:3   t "V(x,y)" with points palette pointsize 3 pointtype 7
#splot "particle10.txt" u 1:2:3   t "n(x,y)" with points palette pointsize 3 pointtype 7


#splot "particle4.txt" u 1:2:3   t "PEPS" with points palette pointsize 3 pointtype 7, 0.2*exp(-(x-2)**2)*exp(-(y-2)**2) t "exp(-x^2)exp(-y^2)" with lines

#{/Symbol @\326\140}2


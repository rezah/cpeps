set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "r_x"   font ",35" textcolor rgb "black"
set ylabel "r_y"  font ",35"  textcolor rgb "black"
set zlabel "<n_i>"  font ",35"  textcolor rgb "black"
set xtics (0,1,2,3,4,5,6);
set ytics (0,1,2,3,4,5,6);
#set ztics (0.0, 0.05,0.1,0.15,0.2,0.25,0.3);
set ztics ( 0.02 ,0.04,0.06,0.08,0.1,0.12);



set tics scale 2

set key box t l 
set key font ",30"
set key spacing 1.5
set key width 2.2


#set output "particle4.eps"
#set output "particle8.eps"
#set output "Potential8.eps"
set output "particle16N10.eps"

#set pm3d
#set dgrid3d 50,50 qnorm 2

set xrange [0:6]
set yrange [0:6]
#set zrange [0:0.12]

set pm3d map
set size square

#splot "particle8.txt" u 1:2:3   t "PEPS" with points palette pointsize 4 pointtype 7, [0:1] 0.02*(sin(x*pi)*sin(2*y*pi)-sin(2*x*pi)*sin(y*pi))**2


#splot [1:0] [0:1] [0:0.12] "particle8.txt" u 1:2:3   t "PEPS" with points palette pointsize 3 pointtype 7, [0:1] [1:0]  0.02*(sin(x*pi)*sin(2*y*pi)-sin(2*x*pi)*sin(y*pi))**2 t "(sin(x)*sin(2y)-sin(2x)*sin(y))^2" with lines



splot "particleN16N10.txt" u 1:2:3   notitle with points palette pointsize 4 pointtype 15

#, 0.09*((sin(x*pi)*sin(2*y*pi)-sin(2*x*pi)*sin(y*pi))**2) t "(sin(x)*sin(2y)-sin(2x)*sin(y))^2" with lines

#{/Symbol @\326\140}2


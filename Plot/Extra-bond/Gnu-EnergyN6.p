set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "1/{/Symbol c}_{o}"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"

set xtics ( "1/4" 0.25,"1/3" 0.333, "1/8" 0.125, "1/16" 0.0625)
#set xtics (0.3,0.2,0.1);
set ytics nomirror (15, 14, 12, 10);
set tics scale 3

set key box t l  
set key font ",30"
set key spacing 2

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

#
G(x)=a+b*(x**(2))+c*(x**(4))
fit   G(x) "N6.txt" u ($1):($2)    via a,b,c
#M(x)=a1+b1*(x**(2))
#fit  [0.03:0.2]  M(x) "N6.txt" u ($1):($2)   via a1,b1#,c
M1(x)=a2+b2*(x**(1))+c2*(x**(3))
fit   M1(x) "N6.txt" u ($1):($2)    via a2,b2,c2



G1(x)=0


set output "ENERGYTLN6.eps"
p [0.01:0.3] [ 10: 14]  10.89507 t "exact" with line lw 4 linetype 7 lc rgb "#5c3566", "N6.txt" u ($1):($2)  t "{/Symbol r}=0.39" with points  pointtype 9 lw 3  ps 5.0 lc rgb "#a40000", x<0.160 ? G(x) : 1/0 notitle lw 4 lc rgb "#f57900" dt 4,x<0.160 ? M1(x) : 1/0 notitle lw 4 lc rgb "#729fcf" dt 4


#x<0.190 ? M(x) : 1/0 notitle lw 4 lc rgb "#729fcf" dt 4

#x<2.0 ? M1(x) : 1/0 notitle lw 4 lc rgb "#75507b" dt 4
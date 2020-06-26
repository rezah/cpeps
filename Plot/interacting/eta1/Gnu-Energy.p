set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "1/N"   font ",30" textcolor rgb "black"
set ylabel ""  font ",30"  textcolor rgb "black"
set ytics nomirror (0.6, 0.4,0.2,0.0, -0.2);

G(x)=a+b*(x**(1))+c*(x**(2))
fit   G(x) "N20.txt" u ($1):($4)    via a,b,c
#M1(x)=a2+b2*(x**(1))+c2*(x**(2))
#fit   M1(x) "N20.txt" u ($1):($4)    via a2,b2,c2


set xtics ( "1/4" 0.25,"1/8" 0.12, "1/16" 0.0625);
set tics scale 3

set key box c r 
set key font ",30"
set key spacing 1.30
set key width 0.90


set output "Density.eps"



p [0:0.3]  [0.61:-0.2] 0.467 t "QMC" with line lw 4 linetype 7 lc rgb "#fa7202","Ne.txt" u 1:2:3  notitle  with errorbars lw 2 , "N20.txt" u ($1):($4)  t  "TN" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#09a5ed" , x<1.230 ? G(x) : 1/0 notitle lw 4 lc rgb "#f57900" dt 4



#pointtype 12 lw 3  ps 4.0 lc rgb "#09a5ed"
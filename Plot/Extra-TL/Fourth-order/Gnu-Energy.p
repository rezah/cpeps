set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "{/Symbol e}"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"

set xtics nomirror (2, 1.5,1,1.25,0.8,0.5,0.25)
#set xtics (0.3,0.2,0.1);
set ytics nomirror (51,48,46, 44,42,40);
set tics scale 3

set key box at 1.23,49
set key font ",25"
set key spacing 1

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(2))+c*(x**(4))
fit   G(x) "N14.txt" u (6/($1+1)):($2)    via a,b,c
#M(x)=a1+b1*(x**(4))
#fit   M(x) "N14.txt" u (6/($1+1)):($2)    via a1,b1#,c
#M1(x)=a2+b2*(x**(2))+c2*(x**(4))
#fit   M1(x) "N14.txt" u (6/($1+1)):($2)    via a2,b2,c2



G1(x)=0


set output "ENERGYTL.eps"
set multiplot

p [0.3:1.3] [40:51] 50.1704 notitle with line lw 3 linetype 7 lc rgb "#006400", "N14.txt" u (6/($1+1)):($2)  t "{/Symbol r}=0.39" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#e3130b", x<1.230 ? G(x) : 1/0 notitle lw 4 lc rgb "#f57900" dt 4

#, x<2.0 ? M1(x) : 1/0 notitle lw 3 lc rgb "#f207a0" dt 3


set origin .17, .12700
set size 0.5,0.50
clear
unset key
unset xlabel
set xtics font "Times-Roman, 25"
set ytics font "Times-Roman, 25"
set xlabel "1/{/Symbol c}_{o}"   font ",25" textcolor rgb "black"
set ylabel ""  font ",25"  textcolor rgb "black"

set xtics ( "1/4" 0.25,"1/3" 0.333, "1/8" 0.125, "1/20" 0.05)
#set xtics (0.3,0.2,0.1);
set ytics nomirror (59, 56, 53, 50, 49,48);
set tics scale 3

set key box at 0.29,51 
set key font ",20"
set key spacing 1.04
set key width 0.5

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

#
G(x)=a+b*(x**(2))+c*(x**(4))
fit [0.03:0.2]  G(x) "N144.txt" u ($1):($2)    via a,b,c
#M(x)=a1+b1*(x**(2))
#fit  [0.03:0.2]  M(x) "N144.txt" u ($1):($2)   via a1,b1#,c
#M1(x)=a2+b2*(x**(2))+c2*(x**(3))
#fit  [0.03:0.2] M1(x) "N144.txt" u ($1):($2)    via a2,b2,c2



G1(x)=0


p [0.01:0.3] [ 48: 55]  49.1598 notitle with line lw 3 linetype 7 lc rgb "#354766", "N144.txt" u ($1):($2)  t "{/Symbol e}=0.67" with points  pointtype 9 lw 3  ps 2.60 lc rgb "#db1846", x<0.190 ? G(x) : 1/0 notitle lw 3 lc rgb "#f57900" dt 4


#,x<0.190 ? M1(x) : 1/0 notitle lw 3 lc rgb "#0fa0bd" dt 3


#x<0.190 ? M(x) : 1/0 notitle lw 4 lc rgb "#729fcf" dt 4

#x<2.0 ? M1(x) : 1/0 notitle lw 4 lc rgb "#75507b" dt 4

unset multiplot



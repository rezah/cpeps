set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 30"
set ytics font "Times-Roman, 30"
set xlabel "{/Symbol e}"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"

set xtics nomirror (2, 1.5,1,1.25,0.8,0.5,0.25)
#set xtics (0.3,0.2,0.1);
set ytics nomirror (51,48,46, 44,42,40);
set tics scale 3

set key box b l  
set key font ",30"
set key spacing 1

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(2))+c*(x**(4))
fit  [0.1:1.1] G(x) "N14.txt" u (6/($1+1)):($2)    via a,b,c
#M(x)=a1+b1*(x**(4))
#fit   M(x) "N14.txt" u (6/($1+1)):($2)    via a1,b1#,c
M1(x)=a2+b2*(x**(2))+c2*(x**(4))
fit   M1(x) "N14.txt" u (6/($1+1)):($2)    via a2,b2,c2

G1(x)=0

set output "ENERGYTL.eps"

p [0.1:1.3] [35:51] 50.1704 notitle with line lw 3 linetype 7 lc rgb "#006400", "N14.txt" u (6/($1+1)):($2)  t "{/Symbol r}=0.39" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#e3130b", x<1.230 ? G(x) : 1/0 notitle lw 4 lc rgb "#f57900" dt 4


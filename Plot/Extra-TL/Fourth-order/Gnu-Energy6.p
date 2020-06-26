set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "{/Symbol e}"   font ",30" textcolor rgb "black"
set ylabel "E"  font ",30"  textcolor rgb "black"

set xtics nomirror (2, 1.5,1,1.25,0.5,0.25)
#set xtics (0.3,0.2,0.1);
set ytics nomirror (10.5,11, 11.2,12);
set tics scale 3

set key box l b  
set key font ",30"
set key spacing 2

g(x,min,max)=( (x>=min && x<=max)? 1.0 : (1/0) )

G(x)=a+b*(x**(2))+c*(x**(4))
fit   G(x) "N6.txt" u (6/($1+1)):($2)    via a,b,c
#M(x)=a1+b1*(x**(4))
#fit   M(x) "N6.txt" u (6/($1+1)):($2)    via a1,b1#,c
M1(x)=a2+b2*(x**(2))+c2*(x**(4))
fit   M1(x) "N6.txt" u (6/($1+1)):($2)    via a2,b2,c2



G1(x)=0


set output "ENERGYTLN6.eps"
p [0.3:1.3] [10.2:11.2] 10.966225 notitle with line lw 4 linetype 7 lc rgb "#006400", "N6.txt" u (6/($1+1)):($2)  t "{/Symbol r}=0.39" with points  pointtype 9 lw 3  ps 5.0 lc rgb "#8B008B", x<2.0 ? G(x) : 1/0 notitle lw 4 lc rgb "#1E90FF" dt 4, x<2.0 ? M1(x) : 1/0 notitle lw 4 lc rgb "#ce5c00" dt 4

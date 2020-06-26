set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "1/{/Symbol c}_{o}"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",40"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics ( "1/4" 0.25,"1/3" 0.333, "1/8" 0.125, "1/16" 0.0625)
set tics scale 4

set key box b r 
set key font ",30"
set key spacing 1.2
set key width 1.20


set output "Density.eps"



p [0.03:0.3] [ -3.50 : -0.1] "FU2.txt" u (1/$1):((log(($2-10.895071280)*(1/10.895071280))/log(10)))  t "{/Symbol r}=0.16" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#f57900","FU3.txt" u (1/$1):((log(($2-26.96744)*(1/26.96744))/log(10)))  t "{/Symbol r}=0.28" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#872071" ,"FU4.txt" u (1/$1):((log(($2-49.1598)*(1/49.1598))/log(10)))  t "{/Symbol r}=0.38" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#204a87" 



#"FU.txt" u (1/$1):((log(($2-9.06647784028)*(1/9.06647784028))/log(10)))  t "{/Symbol r}=0.16, {/Symbol e}=1.0" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ef2929"


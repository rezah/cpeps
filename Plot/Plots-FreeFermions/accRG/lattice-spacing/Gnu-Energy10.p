set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 37"
set ytics font "Times-Roman, 37"
set xlabel "1/{/Symbol c}_{r}"   font ",37" textcolor rgb "black"
set ylabel ""  font ",37"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics ("1/2" 0.5,"1/4" 0.25, "1/6" 0.166,"1/10" 0.1);
set tics scale 4

set key box b r 
set key font ",37"
set key spacing 1.5
#set key width 3.5


set output "RGN10.eps"



p [ 0.01: 0.45] [ -3.3 : -0.80] "EnvAcc44N10.txt" u ($1):((log(($2+10.80687)*(1/10.80687))/log(10)))  t "{/Symbol e}=1.2" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#f57900" ,"EnvAcc66N10.txt" u ($1):((log(($2+15.7627027)*(1/15.7627027))/log(10)))  t "{/Symbol e}=0.86" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#b62d89" ,"EnvAcc88N10.txt" u ($1):((log(($2+17.56658)*(1/17.56658))/log(10)))  t "{/Symbol e}=0.67" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#3465a4" 





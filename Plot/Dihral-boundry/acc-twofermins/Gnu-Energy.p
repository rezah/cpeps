set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "{/Symbol e}"   font ",40" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",40"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics (0.3,0.2,0.1);
set tics scale 4

set key box b r 
set key font ",30"
set key spacing 1.5
#set key width 2.2


set output "Energy.eps"



p [0.01 : 0.25] [ -4 : -0.5] "Energy.txt" u (1/(($1)+1)):((log(($4+69.08723)*(1/69.08723))/log(10)))  t "normal" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#a40000", "Energy.txt" u (1/(($1)+1)):((log(($5+69.08723)*(1/69.08723))/log(10)))  t "antisymmetric" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#204a87"



#"Energy.txt" u (1/(($1)+1)):((log(($2+69.08723)*(1/69.08723))/log(10)))  t "O({{/Symbol e}}^{2}), PEPS" with points  pointtype 5 lw 3  ps 4.0 lc rgb "#f57900", "Energy.txt" u (1/(($1)+1)):((log(($6+69.08723)*(1/69.08723))/log(10)))  t "O({{/Symbol e}}^{2}), exact" with points  pointtype 4 lw 3  ps 4.0 lc rgb "#ef2929",

#"Energy.txt" u (1/(($1)+1)):((log(($3+69.08723)*(1/69.08723))/log(10)))  t "O({{/Symbol e}}^{4}), PEPS" with points  pointtype 6 lw 3  ps 4.0 lc rgb "#ce5c00"

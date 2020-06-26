set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 35"
set ytics font "Times-Roman, 35"
set xlabel "1/D"   font ",35" textcolor rgb "black"
set ylabel "{/Symbol D}E"  font ",35"  textcolor rgb "black"
set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);


#1{/Symbol=30  \264}

set xtics ( "1/4" 0.25,"1/6" 0.1666, "1/8" 0.125, "1/5" 0.2)
set tics scale 4

set key box b r 
set key font ",30"
set key spacing 1.2
set key width 1.20


set output "Density.eps"



p [0.1:0.3] [ -3.80 : -1.50] "FU4.txt" u (1/$1):((log(($2+10.94432)*(1/10.94432))/log(10)))  t "4{/Symbol \264}4" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#f57900","FU6.txt" u (1/$1):((log(($2+26.39124)*(1/26.39124))/log(10)))  t "6{/Symbol \264}6" with points  pointtype 11 lw 3  ps 4.0 lc rgb "#204a87" , "FU8.txt" u (1/$1):((log(($2+48.3264)*(1/48.3264))/log(10)))  t "8{/Symbol \264}8" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#872036" ,"FU10.txt" u (1/$1):((log(($2+76.748)*(1/76.748))/log(10)))  t "10{/Symbol \264}10" with points  pointtype 17 lw 3  ps 4.0 lc rgb "#e810bd"

#, "FU12.txt" u (1/$1):((log(($2+111.65472)*(1/111.65472))/log(10)))  t "12{/Symbol \264}12" with points  pointtype 17 lw 3  ps 4.0 lc rgb "#91072c"  





#"FU.txt" u (1/$1):((log(($2-9.06647784028)*(1/9.06647784028))/log(10)))  t "{/Symbol r}=0.16, {/Symbol e}=1.0" with points  pointtype 13 lw 3  ps 4.0 lc rgb "#ef2929"


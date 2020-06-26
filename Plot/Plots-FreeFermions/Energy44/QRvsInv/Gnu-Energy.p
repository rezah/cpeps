set terminal postscript eps size 5.3,5.3 enhanced color font 'Helvetica,25' lw 3

set xtics font "Times-Roman, 40"
set ytics font "Times-Roman, 40"
set xlabel "iteration"   font ",40" textcolor rgb "black"
set ylabel "E"  font ",40"  textcolor rgb "black"
#set ytics nomirror ( "10^{0}" 0, "10^{-1}" -1, "10^{-2}" -2, "10^{-3}" -3,"10^{-4}" -4, "10^{-5}" -5);



#set xtics (100,200,300);
set tics scale 4

set key box t r 
set key font ",30"
set key spacing 1.5
#set key width 2.5


set output "EnergyD6.eps"

p  "QRD6.txt" u ($1):($3)  t "QR" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#5c3566" ,"InvD6.txt" u (($1)+46):($3)  t "Inv" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#f57900"


#p  "QRD5.txt" u ($1):($3)  t "QR" with points  pointtype 7 lw 3  ps 4.0 lc rgb "#5c3566" ,"InvD5.txt" u (($1)+46):($3)  t "Inv" with points  pointtype 9 lw 3  ps 4.0 lc rgb "#f57900"



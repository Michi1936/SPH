set xrange[-1:11.000000]
set yrange[-1:10.000000]
do for[i=1:t:2]{
print i
plot '1.000000_rigid6_plot.dat' index 0 u 2:3 w p lt 5, '1.000000_rigid6_plot.dat' index i u 2:3 w p lt 7, '1.000000_rigid6_plot.dat' index i+1 u 2:3 w p lt 10
}

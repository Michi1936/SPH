
set xrange[-1:10]	
set yrange[-1:10]
do for[i=1:t:2]{
print i
plot 'plot.dat' index 0 u 2:3 w p lt 5, 'plot.dat' index i u 2:3:8 w p lt 7 lc palette,'plot.dat' index i+1 u 2:3 w p lt 10

}

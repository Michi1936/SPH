set term gif anim delay 2
set out 'rigid6.gif'

set size square
set xrange[-1:30]	
set yrange[-1:30]
do for[i=1:t:2]{
print i
plot 'plot.dat' index 0 u 2:3 w p lt 5, 'plot.dat' index i u 2:3 w p lt 7, 'plot.dat' index i+1 u 2:3 w p lt 8 
}

set out

set term gif anim delay 2
set out 'slope1.gif'

set xrange[-1:11]	
set yrange[-1:21]
do for[i=1:t:1]{
print i
plot 'test13_plot.dat' index 0 u 2:3 w p lt 5, 'test13_plot.dat' index i u 2:3:8 w p lt 7 lc palette
}

set out

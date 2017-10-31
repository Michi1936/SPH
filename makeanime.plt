set term gif anim delay 2
set out 'gif18.gif'

set xrange[-0.5:25]	
set yrange[-0.5:25]
do for[i=0:t:1]{
print i
plot 'plot.dat' index 0 u 2:3 w p, 'plot.dat' index i u 2:3 w p lt 7
}

set out

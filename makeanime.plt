set term gif anim delay 2
set out 'test7_2.gif'

set xrange[-1:17]	
set yrange[-1:21]
do for[i=1:t:1]{
print i
plot 'plot.dat' index 0 u 2:3 w p lt 5, 'plot.dat' index i u 2:3:8 w p lt 7 lc palette
}

set out

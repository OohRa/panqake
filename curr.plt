set title "Example Title"
set title font ",20" norotate

plot for [col=2:5] "currAll.txt" using 1:col with lines
pause -1

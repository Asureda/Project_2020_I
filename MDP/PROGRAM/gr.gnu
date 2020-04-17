set term png
set xlabel"r(r/{/Symbol s})"
set ylabel"g(r)(a.u.)
set key outside

set output 'g_r.png'
plot 'distrib_funct.dat' u 1:2 w l t'g(r)'

#set term eps enhanced
#set output'g_r.eps'
#replot
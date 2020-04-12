set term png
set output 'energy_reduced.png'
set xlabel"t(a.u.)"
set ylabel"energy(a.u.)
set key outside
set style line 1 lc rgb '#E41A1C' lw 1
set style line 2 lc rgb '#440154' lw 1
set style line 3 lc rgb '#21908d' lw 1

plot 'thermodynamics_reduced.dat' u 1:($2) w l ls 1 t'kinetic', '' u 1:($3) w l ls 2 t'pot', '' u 1:($4) w l ls 3 t'tot'

set output 'temperature_reduced.png'
set ylabel"temperature(a.u.)
set style line 4 lc rgb '#e56b5d' lw 1
plot 'thermodynamics_reduced.dat' u 1:5 w l ls 4 t'temperature'

set output 'pressure_reduced.png'
set ylabel"pressure(a.u.)
set style line 5 lc rgb '#0c0887' lw 1
plot 'thermodynamics_reduced.dat' u 1:6 w l ls 5 t'pressure'

#set output 'momentum_reduced.png'
#set ylabel"momentum(a.u.)
#plot 'termodynamics_reduced.dat' u 1:7 w p t'momentum'

#set output 'g_r.png'
#set ylabel"g_r(a.u.)
#plot 'distriv_funct.dat' u 1:2 w l t'g_r'

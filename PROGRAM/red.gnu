set term png
set output 'energy_reduced.png'
set xlabel"t(a.u.)"
set ylabel"energy(a.u.)
set key outside
plot 'thermodynamics_reduced.dat' u 1:($2/256) w p t'kietic','' u 1:($3/256) w p t'pot','' u 1:($4/256) w p t'tot'

set output 'temperature_reduced.png'
set ylabel"temperature(a.u.)
plot 'thermodynamics_reduced.dat' u 1:5 w p t'temperature'

set output 'pressure_reduced.png'
set ylabel"pressure(a.u.)
plot 'thermodynamics_reduced.dat' u 1:6 w p t'pressure'

#set output 'momentum_reduced.png'
#set ylabel"momentum(a.u.)
#plot 'termodynamics_reduced.dat' u 1:7 w p t'momentum'

#set output 'g_r.png'
#set ylabel"g_r(a.u.)
#plot 'distriv_funct.dat' u 1:2 w p t'g_r'
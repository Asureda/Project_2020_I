set term png
set output 'energy_real.png'
set xlabel"t(seconds)"
set ylabel"energy(kJ/mol)"
set key outside
set style line 7 lc rgb '#E41A1C' lw 1 
set style line 8 lc rgb '#440154' lw 1 
set style line 9 lc rgb '#21908d' lw 1 

plot 'thermodynamics_real.dat' u 1:($2) w l ls 7 t'kinetic', '' u 1:($3) w l ls 8 t'pot', '' u 1:($4) w l ls 9 t'tot'

set output 'temperature_real.png'
set ylabel"temperature(Kelvin)"
set style line 10 lc rgb '#e56b5d' lw 1
plot 'thermodynamics_real.dat' u 1:5 w l ls 10 t'temperature'

set output 'pressure_real.png'
set ylabel"pressure(Pascals)"
set style line 11 lc rgb '#0c0887' lw 1
plot 'thermodynamics_real.dat' u 1:6 w l ls 11 t'pressure'

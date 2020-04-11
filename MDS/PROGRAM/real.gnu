set term png
set output 'energy_real.png'
set xlabel"t(seconds)"
set ylabel"energy(kJ/mol)"
set key outside
plot 'thermodynamics_real.dat' u 1:($2) w p t'kietic','' u 1:($3) w p t'pot','' u 1:($4) w p t'tot'

set output 'temperature_real.png'
set ylabel"temperature(Kelvin)"
plot 'thermodynamics_real.dat' u 1:5 w p t'temperature'

set output 'pressure_real.png'
set ylabel"pressure(Pascals)"
plot 'thermodynamics_real.dat' u 1:6 w p t'pressure'

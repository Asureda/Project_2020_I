set term png
set output 'energy_real.png'
set xlabel"t(seconds)"
set ylabel"energy(kJ/mol)"
set key outside
plot 'kin_pot_real.dat' u 1:2 w p t'kietic','kin_pot_real.dat' u 1:3 w p t'pot','total_real.dat' u 1:2 w p t'tot'

set output 'temperature_real.png'
set ylabel"temperature(Kelvin)"
plot 'temp_pres_real.dat' u 1:2 w p t'temperature'

set output 'pressure_real.png'
set ylabel"pressure(Pascals)"
plot 'temp_pres_real.dat' u 1:3 w p t'pressure'
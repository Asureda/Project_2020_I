########################################################
#      GRAPHIC OF THE RADIAL DISTRIBUTION FUNCTION     #
########################################################

# SET TERMINAL
set term png

# RDF (g(r))
set title 'Radial Distribution Function'
set xlabel 'r (Angstrom)'
set ylabel 'g(r)'
set tics font ", 10"
set xrange[1.000000000000000E-002:4.27000000000000]

# Generate output file
set output 'g_r.png'
plot 'distrib_funct.dat' u 1:2 notitle w l lw 1 lc rgb "dark-violet"

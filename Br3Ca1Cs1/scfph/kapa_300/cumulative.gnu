set term png color enh solid
set logscale x
set xlabel 'L (nm)'
set ylabel "Cumulative kappa (W/mK)"
set output "cumulative_300K_kappa.png"
plot "cumulative_300K.dat" using 1:2 w lp

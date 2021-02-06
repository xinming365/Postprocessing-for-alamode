set term png color enh solid
set xlabel "Frequency (cm^{-1})"
set ylabel "{\kappa_L(\omega)}(W/mK/cm^{-1})
set output "spectrum_300K.png"
set xrange [0:350]
plot "pc_rta.kl_spec"  using 2:3 w l lt 2 lw 2

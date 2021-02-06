set term png color enh solid
#set terminal postscript eps enhanced color solid
set xlabel "Temperature (K)"
set ylabel "{/symbol W}_{q}^{2} (cm^{-2})"
set output "Frequency.png"
set key at graph 0.3, graph 0.6
set border lw 2
set mxtics 5
show mxtics
set mytics 5
show mytics
plot "frequency_R.dat"  using 1:($2**2) title "CsCaBr3"  with l lt 2 lw 2 


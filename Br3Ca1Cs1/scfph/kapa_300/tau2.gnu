set term png color enh solid
set size square
set xrange [1:]
set logscale y
set xtics font ',16'  
set ytics font ',16'
set grid lw 2
set border lw 2
#set style arrow 1 head filled size screen 0.025, 30 ,40, ls1
set arrow 1 from 65, 30 to 48, 1.5  head filled size 5, 15, 60 lw 2 lc rgb "skyblue"
set label 1 at 65, 40 center
set label 1 "{/Symbol G} (300 K)" tc lt 3 font 'Times,21' 

set arrow from 160, 20 to 130, 3  head filled size 5, 15, 60 lw 2 lc rgb "light-red"
set label 2 at 160,30 center textcolor rgb "light-red"
set label 2 "{/Symbol G} (600 K)" tc rgb "light-red"  font 'Times,21' 
#set ytics 0.01, 10, 1e3 
set ytics add ("10^{-2}" 0.01) ("10^{-1}" 0.1) ("10^{0}" 1) ("10^{1}" 10) ("10^{2}" 100) ("10^{3}" 1000)
#unset mytics
unset key
set xlabel "Phonon frequency (cm^{-1})" font ',16'
set ylabel "Scattering rates (ps^{-1})" font ',16'
set output "tau_300_600K.png"
plot   "tau300K_10.dat"  using 3:(1/$4) with points pt 6 lc rgb "skyblue" title "300K", "tau600K_10.dat" \
       using 3:(1/$4) with points pt 6 lc rgb "light-red"  title "600K"

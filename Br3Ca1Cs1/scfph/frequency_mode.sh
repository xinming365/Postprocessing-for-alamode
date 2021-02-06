#!/bin/bash
M=0.588792
R=1.515039
grep "${R}" pc_scfph.scph_bands | awk -F "    " '{print $2}' > temp.dat
grep "1.515039" pc_scfph.scph_bands | awk -F "    " '{print $3}'| awk '{print $2}' > frequency.dat
paste temp.dat frequency.dat > frequency_R.dat
rm temp.dat frequency.dat
 

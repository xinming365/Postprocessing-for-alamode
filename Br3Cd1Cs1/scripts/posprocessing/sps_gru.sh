#!/bin/bash
rm sps.dat gru.dat
T=300
mode=15
q_num=1728 # 12*12*12=1728

grep "   ${T}  " pc_dos.sps_Bose | awk -v val=$j '{print val, $3,$5,$6}' > sps.dat
for ((i=0;i<${q_num};i=i+1))
do
grep -A ${mode} "knum =" pc_dos.gru_all | head -$(ps -ef | echo "(${mode}+1)*$i"|bc) | tail -15 | awk -v val=$j '{print val, $3,$4}' >> gru.dat
done


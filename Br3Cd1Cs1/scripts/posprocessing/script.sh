rm kappa.dat
for ((i=100;i<=900;i=i+100)) 
do 
grep -A 2 Temperature kapa_${i}/Rb3AuO_rta.kl | tail -1 >> kappa.dat
done


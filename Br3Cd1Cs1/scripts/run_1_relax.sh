#!/bin/bash
#SBATCH -p v3_64
#SBATCH -N 2
#SBATCH -n 40
source /public1/home/sc30830/software/alamode/alamode-1.1.0/env.sh
source /public1/soft/other/vasp/cn-module-vasp.5.4.4.sh
module  unload intel/17.0.7
module load  mpi/intel/5.0.3.049
module load vasp/intel-17/vasp544
module load anaconda/3-phonopy-phono3py

# This script is the version for slurm system.

#VASP="/home/nijun/vasp.5.3.2/vasp.5.3/vasp.nosoc"
#vasppot="/home/nijun/apps/vasp/vasppot/vasppot_paw_pbe.52.sh"
vasppot="/public1/home/sc30830/pseudo/vasppot/vasppot_paw_pbe.54.sh"

#-----------------------------------------------------------------------#
element_num=3   # No. of chemical elements in your system !????????????????????????
echo ${element_num} > element_num.dat
#-----------------------------------------------------------------------#
mkdir relax2
cd ./relax2
cat > INCAR << eof
SYSTEM=Br3CsCd
ISTART=0
ICHARG=2
ENCUT=650

ISMEAR=0
SIGMA=0.01
EDIFF=1E-8
NELMIN=5
EDIFFG=-1E-6
ISIF=3
IBRION=2
POTIM=0.1
NSW=500
#ISPIN=2
PREC=Accurate

GGA=PS
IALGO=38
LWAVE=.FALSE.
LCHARG=.FALSE.
LREAL=.FALSE.
ADDGRID=.TRUE.
#NPAR=4
eof
cat > KPOINTS << eof
A
0
G
12 12 12
 0  0  0
eof
cat > POSCAR <<eof
Br3Cd1Cs1   [CUB,CUB,cP5] (STD_PRIM doi:
1.000000
   5.60081618344005   0.00000000000000   0.00000000000000
   0.00000000000000   5.60081618344005   0.00000000000000
  -0.00000000000000   0.00000000000000   5.60081618344005
Br Cd Cs
3 1 1 
Direct(5) [A3B1C1] 
   0.00000000000000   0.50000000000000   0.50000000000000      
   0.50000000000000   0.00000000000000   0.50000000000000      
   0.50000000000000   0.50000000000000   0.00000000000000      
   0.50000000000000   0.50000000000000   0.50000000000000      
   0.00000000000000   0.00000000000000   0.00000000000000  
eof
#vasppot_paw61_pbe.sh Li_sv H
${vasppot} Br  Cd  Cs_sv
#${MPIVASP} > relax.out
srun vasp_std > relax.out

#condi=$(cat relax.out | tail -1 | awk '{print $1}')
#while [ "${condi}" != "reached" ]
#do 
#cp CONTCAR POSCAR
#i${MPIUSE} ${VASP5_nosoc} > relax.out
#${MPIVASP} > relax.out
#condi=$(cat relax.out | tail -1 | awk '{print $1}')
#done

cp CONTCAR POSCAR
cat > INCAR << eof
SYSTEM=Br3CaCs
ISTART=0
ICHARG=2
ENCUT=650

ISMEAR=0
SIGMA=0.01
EDIFF=1E-6
NELMIN=5
#LCALCEPS=.TRUE.  # Static dielectric tensor, Born effective charges, and piezoelectric constants for DFT, HSE06, etc.
LEPSILON=.TRUE.  # Static dielectric tensor, Born effective charges, and piezoelectric constants for standard DFT.
LPEAD=.TRUE.     # Derivative of cell-periodic part of wave function w.r.t. Bloch vector, use with LEPSILON=.TRUE. or LOPTICS=.TRUE.
IBRION=8
#ISPIN=2
PREC=Accurate

GGA=PS
IALGO=38
LREAL=.FALSE.
LWAVE=.FALSE.
LCHARG=.FALSE.
ADDGRID=.TRUE.
#NPAR=4
NELM=100
eof
cat > KPOINTS << eof
A
0
G
12 12 12
 0  0  0
eof
#${MPIUSE} ${VASP_nosoc} > scf_born.out
srun vasp_std > scf_born.out

grep -A 4 "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)" OUTCAR | tail -3 |awk -v val=$j '{print val,$1,$2,$3}' > BORN

ion_1=$(grep -B 1 Direct POSCAR | head -1 | awk '{print $1}')
ion_2=$(grep -B 1 Direct POSCAR | head -1 | awk '{print $2}')
ion_3=$(grep -B 1 Direct POSCAR | head -1 | awk '{print $3}')
ion_4=$(grep -B 1 Direct POSCAR | head -1 | awk '{print $4}')
if [ "$element_num" = "1" ];then
ion=${ion_1}
elif [ "$element_num" = "2" ];then
ion=$(echo "${ion_1}+${ion_2}"|bc)
elif [ "$element_num" = "3" ];then
ion=$(echo "${ion_1}+${ion_2}+${ion_3}"|bc)
elif [ "$element_num" = "4" ];then
ion=$(echo "${ion_1}+${ion_2}+${ion_3}+${ion_4}"|bc)
else
echo "The No. of chemical elements is wrong."
fi

num=$(echo "$ion*4+1" | bc)
for ((i=1;i<${num};i=i+4))
do
grep -A $(echo "$i+4"|bc) "BORN EFFECTIVE CHARGES (in e, cummulative output)" OUTCAR | tail -3 |awk -v val=$j '{print val,$2,$3,$4}' >> BORN
#grep -A $(ps -ef |echo "$i+4"|bc) "BORN EFFECTIVE CHARGES (including local field effects)" OUTCAR | tail -3 |awk -v val=$j '{print val,$2,$3,$4}' >> BORN
done

cp POSCAR POSCAR_pc  
phonopy --symmetry POSCAR   # produce unite cell
cp BPOSCAR POSCAR           # use unite cell to 
phonopy -d --dim="2 2 2"    # generate supercell  !????????????????????????


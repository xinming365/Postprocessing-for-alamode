#!/bin/bash
#SBATCH -p v3_64
#SBATCH -N 2
#SBATCH -n 40
source /public1/soft/modules/module.sh
source /public1/home/sc30830/software/alamode/alamode-1.1.0/env.sh
source /public1/soft/other/vasp/cn-module-vasp.5.4.4.sh
module  unload intel/17.0.7
module unload anaconda
module load python/3.6.5
module load  mpi/intel/5.0.3.049
module load vasp/intel-17/vasp544
module load  anaconda/3-phonopy-phono3py

#VASP="/home/nijun/apps/vasp/vasp.5.4.1/bin/vasp_std"
vasppot="/public1/home/sc30830/pseudo/vasppot/vasppot_paw_pbe.54.sh"


#find ./md -name "POSCAR_*" -print > filename
num=$(awk 'END{print NR}' filename)
#
#for (( i=1; i <= ${num}; i=i+1 ))
#do 
#mkdir job-$i
#mv ./md/POSCAR_$i job-$i/POSCAR
#cd ./job-$i
#cat > INCAR << eof
#SYSTEM=Br3CaCs
#ISTART=0
#ICHARG=2
#ENCUT=650
#
#ISMEAR=0
#SIGMA=0.01
#EDIFF=1E-8
#NELMIN=5
#IBRION=-1
##ISPIN=2
#PREC=Accurate
#
#GGA=PS
#IALGO=38
#LREAL=.FALSE.
#LWAVE=.FALSE.
#LCHARG=.FALSE.
#ADDGRID=.TRUE.
#NPAR=4
##SYMPREC=0.001
#NELM=100
#eof
#cat > KPOINTS << eof
#A
#0
#G
#2 2 2
#0 0 0
#eof
##vasppot_paw61_pbe.sh Li_sv H
#${vasppot} Br Ca_sv Cs_sv
##${MPIUSE} ${VASP_nosoc} > scf.out
#srun vasp_std > scf.out
#cd ..
#done
#
rm disp-force-an.dat
for (( i=1; i <= ${num}; i=i+1 ))
do
python3 ~/software/alamode/alamode-1.1.0/tools/extract.py --VASP=POSCAR job-$i/vasprun.xml >> disp-force-an.dat
#python ~/apps/alamode/alamode-v.1.1.0/tools/extract.py --VASP=POSCAR job-$i/vasprun.xml >> disp-force-an.dat
done

#-----------------------------------------------------------------------
element_num=$(cat element_num.dat)
element=$(grep -B 2 Direct POSCAR | head -1)
ion=$(cat ion.dat)
## basis vestor of unite cell ##
A_1=$(grep -B 5 Direct POSCAR | head -1 | tail -1)
A_2=$(grep -B 5 Direct POSCAR | head -2 | tail -1)
A_3=$(grep -B 5 Direct POSCAR | head -3 | tail -1)
## basis vestor of unite cell ##
#------------------------------------------------------------------------

## create alm.in with MODE=optimize and NORDER=5 ##
## for cross-validation for obtaiining the suitable L1_ALPHA value ##
cat > alm.in << eof
&general
  PREFIX = uc 
  MODE = optimize         # <-- here
  NAT = ${ion}; NKD = ${element_num}
  KD = ${element}
  PERIODIC = 1 1 1
/
&interaction
  NORDER = 5               # 1: harmonic, 2: cubic, ..
  NBODY = 2 3 3 2 2
/
&optimize
  LMODEL= enet
  NDATA = ${num}
  DFSET = disp-force-an.dat
  FC2XML = uc_2nd.xml  # use the 2nd IFC produced by run_2*.sh step
  L1_RATIO = 1.0            # LASS . beta value in 'Elastic-net regression' section of manual
  CV = 4
  CV_MINALPHA = 1e-10
  CV_MAXALPHA = 0.02
  CV_NALPHA = 100
/
&cutoff
  *-*  None None  18.5 15.1  11.0  # use 3rd, 2nd, 1st nearest neighbor shell for 4, 5, 6 order terms !?????????????????????
  Br-Cd None None 15.6  11.7  5.3
  Br-Cs None None 18.0  12.8  7.4
  Cd-Cs None None 18.0  12.8  9.0
  Br-Br None None 12.8  10.4  7.4
  Cd-Cd None None 18.0  14.7  10.4
  Cs-Cs None None 18.0  14.7  10.4 
/
&cell
  1.8897261254578282       # factor in Bohr unit
  ${A_1}
  ${A_2}
  ${A_3}
/
&position
eof
cat posi.in >> alm.in
alm alm.in > alm3.log

## create alm.in with MODE=optimize and NORDER=5 ##
## for cross-validation for obtaiining the suitable L1_ALPHA value ##
cat > alm.in << eof
&general
  PREFIX = uc
  MODE = optimize          # <-- here
  NAT = ${ion}; NKD = ${element_num}
  KD = ${element}
  PERIODIC = 1 1 1
/
&interaction
  NORDER = 5               # 1: harmonic, 2: cubic, ..
  NBODY = 2 3 3 2 2
/
&optimize
  LMODEL= enet
  NDATA = ${num}
  DFSET = disp-force-an.dat
  FC2XML = uc_2nd.xml      # use the 2nd IFC produced by run_2*.sh step
  L1_RATIO = 1.0           # LASS . beta value in 'Elastic-net regression' section of manual
  L1_ALPHA = 1.41421e-06        # alpha value in 'Elastic-net regression' section of manual     
  CV = 0
/
&cutoff
  *-*  None None  18.5 15.1  11.0  # use 3rd, 2nd, 1st nearest neighbor shell for 4, 5, 6 order terms !?????????????????????
  Br-Cd None None 15.6  11.7  5.3
  Br-Cs None None 18.0  12.8  7.4
  Cd-Cs None None 18.0  12.8  9.0
  Br-Br None None 12.8  10.4  7.4
  Cd-Cd None None 18.0  14.7  10.4
  Cs-Cs None None 18.0  14.7  10.4   
/
&cell
  1.8897261254578282       # factor in Bohr unit
  ${A_1}
  ${A_2}
  ${A_3}
/
&position
eof
cat posi.in >> alm.in
alm alm.in > alm4.log

########################### IMPORTANT !!! ############################
### At the end, you should rerun the last step use the recomended  ###
### value of L1_ALPHA produced by the cross-validation step        ###
######################################################################


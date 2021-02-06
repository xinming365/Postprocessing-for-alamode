#!/bin/bash
#SBATCH -p v3_64
#SBATCH -N 1
#SBATCH -n 24
source /public1/soft/modules/module.sh
source /public1/home/sc30830/software/alamode/alamode-1.1.0/env.sh
source /public1/soft/other/vasp/cn-module-vasp.5.4.4.sh
module  unload intel/17.0.7
module load python/3.6.5
module load  mpi/intel/5.0.3.049
module load vasp/intel-17/vasp544
module load  anaconda/3-phonopy-phono3py

#VASP="/home/nijun/apps/vasp/vasp.5.4.1/bin/vasp_std"
vasppot="/public1/home/sc30830/pseudo/vasppot/vasppot_paw_pbe.54.sh"


find ./md -name "POSCAR_*" -print > filename
num=$(awk 'END{print NR}' filename)

for (( i=1; i <= ${num}; i=i+1 ))
do 
mkdir job-$i
mv ./md/POSCAR_$i job-$i/POSCAR
cd ./job-$i
cat > INCAR << eof
SYSTEM=Br3CdCs
ISTART=0
ICHARG=2
ENCUT=650

ISMEAR=0
SIGMA=0.01
EDIFF=1E-8
NELMIN=5
IBRION=-1
#ISPIN=2
PREC=Accurate

GGA=PS
IALGO=38
LREAL=.FALSE.
LWAVE=.FALSE.
LCHARG=.FALSE.
ADDGRID=.TRUE.
NPAR=4
#SYMPREC=0.001
NELM=100
eof
cat > KPOINTS << eof
A
0
G
2 2 2
0 0 0
eof
#vasppot_paw61_pbe.sh Li_sv H
${vasppot} Br  Cd  Cs_sv
#${MPIUSE} ${VASP_nosoc} > scf.out
cat > sub.sh << eof
#!/bin/bash
#SBATCH -p v3_64
#SBATCH -N 1
#SBATCH -n 24
source /public1/soft/other/vasp/cn-module-vasp.5.4.4.sh
module  unload intel/17.0.7
module load  mpi/intel/5.0.3.049
module load vasp/intel-17/vasp544
srun vasp_std > scf.out
eof
sbatch sub.sh
cd ..
done

#rm disp-force-an.dat
#for (( i=1; i <= ${num}; i=i+1 ))
#do
#python3 ~/software/alamode/alamode-1.1.0/tools/extract.py --VASP=POSCAR job-$i/vasprun.xml >> disp-force-an.dat
##python ~/apps/alamode/alamode-v.1.1.0/tools/extract.py --VASP=POSCAR job-$i/vasprun.xml >> disp-force-an.dat
#done
#
##-----------------------------------------------------------------------
#element_num=$(cat element_num.dat)
#element=$(grep -B 2 Direct POSCAR | head -1)
#ion=$(cat ion.dat)
### basis vestor of unite cell ##
#a_1=$(grep -B 5 Direct POSCAR | head -1 | tail -1)
#a_2=$(grep -B 5 Direct POSCAR | head -2 | tail -1)
#a_3=$(grep -B 5 Direct POSCAR | head -3 | tail -1)
### basis vestor of unite cell ##
#------------------------------------------------------------------------
#
## create alm.in with MODE=optimize and NORDER=5 ##
## for cross-validation for obtaiining the suitable L1_ALPHA value ##
#cat > alm.in << eof
#&general
#  prefix = uc 
#  mode = optimize         # <-- here
#  nat = ${ion}; NKD = ${element_num}
#  kd = ${element}
#  periodic = 1 1 1
#/
#&interaction
#  norder = 5               # 1: harmonic, 2: cubic, ..
#  nbody = 2 3 3 2 2
#/
#&optimize
#  lmodel= enet
#  ndata = ${num}
#  dfset = disp-force-an.dat
#  fc2xml = uc_2nd.xml  # use the 2nd IFC produced by run_2*.sh step
#  l1_ratio = 1.0            # LASS . beta value in 'Elastic-net regression' section of manual
#  cv = 4
#  cv_minalpHA = 1e-10
#  cv_maxalpHA = 0.02
#  cv_nalpha = 100
#/
#&cutoff
#  *-*  none None 9.5 5.4 3.8  # use 6th, 2nd, 1st nearest neighbor shell for 4, 5, 6 order terms !?????????????????????
#/
#&cell
#  1.8897261254578282       # factor in Bohr unit
#  ${a_1}
#  ${a_2}
#  ${a_3}
#/
#&position
#eof
#cat posi.in >> alm.in
#${alm} alm.in > alm3.log
#
### create alm.in with MODE=optimize and NORDER=5 ##
### for cross-validation for obtaiining the suitable L1_ALPHA value ##
#cat > alm.in << eof
#&general
#  prefix = uc
#  mode = optimize          # <-- here
#  nat = ${ion}; NKD = ${element_num}
#  KD = ${element}
#  PERIODIC = 1 1 1
#/
#&interaction
#  NORDER = 5               # 1: harmonic, 2: cubic, ..
#  NBODY = 2 3 3 2 2
#/
#&optimize
#  LMODEL= enet
#  NDATA = ${num}
#  DFSET = disp-force-an.dat
#  FC2XML = uc_2nd.xml      # use the 2nd IFC produced by run_2*.sh step
#  L1_RATIO = 1.0           # LASS . beta value in 'Elastic-net regression' section of manual
#  L1_ALPHA = 2.09128e-07        # alpha value in 'Elastic-net regression' section of manual     
#  CV = 0
#/
#&cutoff
#  *-*  None None 9.5 5.4 3.8  # use 6th, 2nd, 1st nearest neighbor shell for 4, 5, 6 order terms !?????????????????????
#/
#&cell
#  1.8897261254578282       # factor in Bohr unit
#  ${A_1}
#  ${A_2}
#  ${A_3}
#/
#&position
#eof
#cat posi.in >> alm.in
#${alm} alm.in > alm4.log

########################### IMPORTANT !!! ############################
### At the end, you should rerun the last step use the recomended  ###
### value of L1_ALPHA produced by the cross-validation step        ###
######################################################################


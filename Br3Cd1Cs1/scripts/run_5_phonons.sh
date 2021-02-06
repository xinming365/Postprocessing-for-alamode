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


#!/bin/bash
date > output


#-----------------------------------------------------------
element_num=$(cat element_num.dat)
element=$(grep -B 2 Direct POSCAR | head -1)
mass=(79.904 112.41 132.91) # the mass of Br Cd Cs

## basis vestor of primitive cell ##
a_1=$(grep -B 5 Direct relax/POSCAR_pc | head -1 | tail -1)
a_2=$(grep -B 5 Direct relax/POSCAR_pc | head -2 | tail -1)
a_3=$(grep -B 5 Direct relax/POSCAR_pc | head -3 | tail -1)
## basis vestor of primitive cell ##
#-----------------------------------------------------------

mkdir phonons
cp uc_2nd.xml ./phonons
cp uc.xml ./phonons
cp ./relax/BORN ./phonons
cd ./phonons

### phonon dispersion ####
cat > anphon.in << eof
&general
  PREFIX = pc_disp
  MODE = phonons
  NKD = ${element_num}
  KD = ${element}
  MASS = ${mass[@]}      # ????????????????????????????
  FCSXML = uc.xml
  FC2XML = uc_2nd.xml
  NONANALYTIC = 3
  #NA_SIGMA = 0.3
  BORNINFO = BORN
  BORNSYM = 1
  BCONNECT = 2
  EMIN = 0
  EMAX = 1200
  DELTA_E = 2
/
&cell
  1.8897261254578282 # factor in Bohr unit
   ${a_1}
   ${a_2}
   ${a_3}
/
&kpoint
  1  # KPMODE = 1: line mode                # ????????????????????????????????????????
 G 0.0 0.0 0.0  X 0.0 0.5 0.0  50
 X 0.0 0.5 0.0  M 0.5 0.5 0.0  50
 M 0.5 0.5 0.0  G 0.0 0.0 0.0  50 
 G 0.0 0.0 0.0  R 0.5 0.5 0.5  50
 R 0.5 0.5 0.5  M 0.0 0.5 0.0  50
/
&analysis
  GRUNEISEN = 1
  PRINTEVEC = 1
  PRINTXSF = 1
  PRINTVEL = 1
  PRINTPR = 1
/
eof
#${MPIUSE} ${anphon} anphon.in > anphon_disp.log
anphon anphon.in > anphon_disp.log

### phonon dos pdos ######
cat > anphon.in << eof
&general
  PREFIX = pc_dos
  MODE = phonons
  NKD = ${element_num}
  KD = ${element}
  MASS = ${mass[*]}     # ????????????????????????????
  FCSXML = uc.xml
  FC2XML = uc_2nd.xml
  NONANALYTIC = 3
  #NA_SIGMA = 0.3
  BORNINFO = BORN
  BORNSYM = 1
  EMIN = 0
  EMAX = 1200
  DELTA_E = 2
  ISMEAR = -1
  EPSILON = 2
  BCONNECT = 2
/
&cell
  1.8897261254578282 # factor in Bohr unit
   ${a_1}
   ${a_2}
   ${a_3}
/
&kpoint
  2  # KPMODE = 2: Uniform k grid for phonon DOS and thermal conductivity
  12 12 12
/
&analysis
  GRUNEISEN = 1
  PRINTVEL = 1
  PRINTMSD = 1
  PDOS = 1
  TDOS = 1
  SPS = 2
  ANIME = 0 0 0
  ANIME_CELLSIZE = 2 2 2
  ANIME_FORMAT = xsf
  #FE_BUBBLE = 1
/
eof
#${MPIUSE} ${anphon} anphon.in > anphon_dos.log
anphon anphon.in > anphon_dos.log

### thermal conductivite by RTA ###
cat > anphon.in << eof
&general
  PREFIX = pc_rta
  MODE = RTA
  NKD = ${element_num}
  KD = ${element}
  MASS = ${mass[*]}      # ????????????????????????????
  FCSXML = uc.xml
  FC2XML = uc_2nd.xml
  NONANALYTIC = 3
  #NA_SIGMA = 0.3
  BORNINFO = BORN
  BORNSYM = 1
  EMIN = 0
  EMAX = 1200
  DELTA_E = 2
  TMIN = 100
  TMAX = 700
  DT = 50
  ISMEAR = -1
  EPSILON = 2
  BCONNECT = 2
  TRISYM = 1
/
&cell
  1.8897261254578282 # factor in Bohr unit
   ${a_1}
   ${a_2}
   ${a_3}
/
&kpoint
  2  # KPMODE = 2: Uniform k grid for phonon DOS and thermal conductivity
  12 12 12
/
&analysis
  KAPPA_SPEC = 1
  ISOTOPE = 2
/
eof
#${MPIUSE} ${anphon} anphon.in > anphon_rta.log
anphon anphon.in > anphon_rta.log

#~/apps/alamode/alamode-v.1.1.0/tools/analyze_phonons.py --calc kappa_boundary --size 1.0e+6 pc_rta.result > pc_rta_boundary_1mm.kl
python3 ~/software/alamode/alamode-1.1.0/tools/analyze_phonons.py --calc kappa_boundary --size 1.0e+6 pc_rta.result > pc_rta_boundary_1mm.kl
#############################----END----#####################################

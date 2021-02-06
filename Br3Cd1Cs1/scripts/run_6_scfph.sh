#!/bin/bash
#SBATCH -p v3_64
#SBATCH -N 2
#SBATCH -n 40
source /public1/soft/modules/module.sh
source /public1/home/sc30830/software/alamode/alamode-1.1.0/env.sh
source /public1/soft/other/vasp/cn-module-vasp.5.4.4.sh
module  unload intel/17.0.7
module load python/3.6.5
module load  mpi/intel/5.0.3.049
module load vasp/intel-17/vasp544
module load  anaconda/3-phonopy-phono3py
module load gcc/5.4.0

#VASP="/home/nijun/apps/vasp/vasp.5.4.1/bin/vasp_std"
vasppot="/public1/home/sc30830/pseudo/vasppot/vasppot_paw_pbe.54.sh"
dfc2=~/software/alamode/alamode-1.1.0/tools/dfc2
#----------------------------------------------------------
element_num=$(cat element_num.dat)
element=$(grep -B 2 Direct POSCAR | head -1)
mass=(79.904 112.41 132.91) # the mass of Br Cd Cs

## basis vestor of primitive cell ##
a_1=$(grep -B 5 Direct relax/POSCAR_pc | head -1 | tail -1)
a_2=$(grep -B 5 Direct relax/POSCAR_pc | head -2 | tail -1)
a_3=$(grep -B 5 Direct relax/POSCAR_pc | head -3 | tail -1)
## basis vestor of primitive cell ##
#----------------------------------------------------------

mkdir scfph2
cp uc_2nd.xml ./scfph2
cp uc.xml ./scfph2
cp ./relax/BORN ./scfph2
cd ./scfph2

### scf-phonon dispersion ####
cat > anphon.in << eof
&general
  PREFIX = pc_scfph
  MODE = SCPH
  NKD = ${element_num}
  KD = ${element}
  MASS = ${mass[*]}    
  FCSXML = uc.xml
  FC2XML = uc_2nd.xml
  NONANALYTIC = 2
  #NA_SIGMA = 0.3
  BORNINFO = BORN
  BORNSYM = 1
  TMIN = 100
  TMAX = 700
  DT = 50
  EMIN = 0
  EMAX = 1200
  DELTA_E = 2
  BCONNECT = 2
/
&scph
  KMESH_INTERPOLATE = 2 2 2  # same as supercell size, e.g. for 4*4*4 supercell, use a KMESH_INTERPOLATE of 4 4 4.
  KMESH_SCPH = 4 4 4         # equal to or a multiple of the number of KMESH_INTERPOLATE in the same direction.
  SELF_OFFDIAG = 1
  MIXALPHA = 0.1
  MAXITER = 1000
  RESTART_SCPH = 0
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
eof
#${MPIUSE} ${anphon} anphon.in > anphon_scfph_band.log
anphon anphon.in > anphon_scfph_band.log

### scf-phonon dos ####
cat > anphon.in << eof
&general
  PREFIX = pc_scfph
  MODE = SCPH
  NKD = ${element_num}
  KD = ${element}
  MASS = ${mass[*]}  
  FCSXML = uc.xml
  FC2XML = uc_2nd.xml 
  NONANALYTIC = 2
  #NA_SIGMA = 0.3
  BORNINFO = BORN
  BORNSYM = 1
  TMIN = 100
  TMAX = 700
  DT = 50
  EMIN = 0
  EMAX = 1200
  DELTA_E = 2
  BCONNECT = 2
/
&scph
  KMESH_INTERPOLATE = 2 2 2  # same as supercell size, e.g. for 4*4*4 supercell, use a KMESH_INTERPOLATE of 4 4 4.
  KMESH_SCPH = 4 4 4         # equal to or a multiple of the number of KMESH_INTERPOLATE in the same direction.
  SELF_OFFDIAG = 1
  MIXALPHA = 0.1
  MAXITER = 1000
  RESTART_SCPH = 1
/
&cell
  1.8897261254578282 # factor in Bohr unit
   ${a_1}
   ${a_2}
   ${a_3}
/
&kpoint
  2  # KPMODE = 2 unifor k grid for phonon DOS and thermal conductivity
  12 12 12
/
&analysis
  DOS = 1
  PRINTMSD = 1
/
eof
#${MPIUSE} ${anphon} anphon.in > anphon_scfph_dos.log
anphon anphon.in > anphon_scfph_dos.log

############### obtain the an-IFCs corrected IFC2 ## 
for ((i = 100; i <= 700; i = i + 50))
do
${dfc2} <<TextForInput
  uc_2nd.xml
  uc_2nd_${i}K.xml
  pc_scfph.scph_dfc2
  ${i}
TextForInput
done

### use an-IFCs and the corrected IFC2 to recalculate
### phonon dispersion, DOS, and rta-Kl 

for ((i = 100; i <= 700; i = i + 50))
do
mkdir kapa_${i}
cp uc_2nd_${i}K.xml ./kapa_${i}
cp uc.xml ./kapa_${i}
cp BORN ./kapa_${i}
cd ./kapa_${i}
### phonon dispersion ####
cat > anphon.in << eof
&general
  PREFIX = pc_disp    
  MODE = phonons
  NKD = ${element_num}
  KD = ${element}
  MASS = ${mass[*]}
  FCSXML = uc.xml
  FC2XML = uc_2nd_${i}K.xml
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
  MASS = ${mass[*]}
  FCSXML = uc.xml
  FC2XML = uc_2nd_${i}K.xml
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
  MASS = ${mass[*]}
  FCSXML = uc.xml
  FC2XML = uc_2nd_${i}K.xml
  NONANALYTIC = 3
  #NA_SIGMA = 0.3
  BORNINFO = BORN
  BORNSYM = 1
  EMIN = 0
  EMAX = 1200
  DELTA_E = 2
  TMIN = ${i}
  TMAX = ${i}
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
cd ..
done


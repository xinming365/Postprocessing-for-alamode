#!/bin/bash
#SBATCH -p v3_64
#SBATCH -N 1
#SBATCH -n 24
source /public1/home/sc30830/software/alamode/alamode-1.1.0/env.sh
source /public1/soft/other/vasp/cn-module-vasp.5.4.4.sh
module  unload intel/17.0.7
module load python/3.6.5
module load  mpi/intel/5.0.3.049
module load vasp/intel-17/vasp544

vasppot="/public1/home/sc30830/pseudo/vasppot/vasppot_paw_pbe.54.sh"


# do MD & obtain the supercell per 50 steps with 2 fs time step #
mkdir md
cp POSCAR ./md
cp ion.dat ./md
cd ./md
cat > INCAR << eof
SYSTEM=Br3CdCs  #????????????????????????????????
#Genieral:
ISTART = 0
ICHARG = 2
ISMEAR = 0      ! -1 is also possible
SIGMA  = 0.1    ! or 0.05
GGA=PS

#Precision:
PREC   = Normal
EDIFF  = 1E-4

#Electronic Step:
NELM   = 300
NELMIN = 6      ! between 6 and 8.
ENCUT= 650
ALGO=Fast       ! Use ALGO=Very Fast for large molecular dynamics runs.
ISYM=0
MAXMIX=40       ! about three times of the number of iterations in the first iteration.
NPAR=4

#MD 
IBRION = 0      ! do MD 
NSW =  3000     
NWRITE = 0 
TEBEG = 300
TEEND = 300         
SMASS = 2.0     ! -3 is micro-canonical ensemble. -1 do annealing . > 0 NVT .
                ! Usually use 2 for NVT, increase it will decrease the Nose-frequency (in Hz)
                ! The Nose-frequency in OUTCAR should be similar to phonon frequency.
NBLOCK = 50     ! NBLOCK steps write one time XDATCAR.   
POTIM = 2       ! timestep 1 fs 1.0-3.0
IWAVPR=11       ! to resolve warning of "Information: wavefunction orthogonal band ** "

#Output
LWAVE=.FALSE.
LCHARG=.FALSE.
LREAL=Auto
eof
cat > KPOINTS << eof
A
0
G
2  2  2 
0  0  0
eof
#vasppot_paw61_pbe.sh Li_sv H
${vasppot} Br Cd  Cs_sv   # !??????????????????????????
#${MPIUSE} ${VASP_nosoc} > md.out
srun vasp_std > md.out

num=60   ## there are 60 supercells in XDATCAR file
ion=$(cat ion.dat)  ## there are 40 atoms in each supsercell
lattice=$(cat POSCAR | head -3 | tail -1 | awk '{printf("%.5f\n",$1)}')  ## lattice constane of the supercell

## change 3 places labeled by '!#' in random.f90 file for differernt supercell system ##
cat > random.f90 << eof
Program suijishu
  implicit none
  integer counter
  Real(kind=8):: x(3),y(3),radious,lattice,m1,m2,m3
  open(unit=1,File='random.dat')
  open(unit=2,File='direct.dat')
  counter=1
  lattice=${lattice}
  call random_seed()
  do while (counter.LE.${num}*${ion})      
    call random_number(x)
    x=x*0.4-0.2
    radious=dsqrt(x(1)**2+x(2)**2+x(3)**2)
    if (dabs(radious-0.1)<=1e-6) then
      write(1,'(4F15.8)') x,radious
      m1=x(1)/lattice  !# For SC supercell, change the Car
      m2=x(2)/lattice  !# coordinate to Direct one! For other
      m3=x(3)/lattice  !# supercell, these lines should be changed
      write(2,'(4F15.8)')m1,m2,m3
      counter=counter+1
      write(*,*)counter-1
    end if
  end do
  close(1)
End Program
eof
ifort random.f90
./a.out

cat > produce_POSCAR.f90 << eof
Program POSCAR_all
  implicit none
  character*50::fname,filename,line_1,line_2
  integer :: i,j
  Real(kind=8):: lattice_scale,a1(3),a2(3),a3(3),x(3),d(3)
  open(unit=1,File='XDATCAR')
  open(unit=2,File='direct.dat')

  read(1,*)
  read(1,*)lattice_scale
  read(1,*)a1
  read(1,*)a2
  read(1,*)a3
  read(1,100)line_1  
  100 format(A25)
  read(1,200)line_2
  200 format(A23)
 
  do i=1,${num}           
    write(fname,'(i4)')i
    filename="POSCAR_"//trim(adjustL(fname))
    open(unit=10,file=filename)
    write(10,*)"supercell",i
    write(10,*)lattice_scale
    write(10,'(3F15.8)')a1
    write(10,'(3F15.8)')a2
    write(10,'(3F15.8)')a3
    write(10,*)line_1
    write(10,*)line_2
    write(10,'(6A)')"Direct"
    read(1,*)
    do j=1,${ion}
      read(1,*)x  
      read(2,*)d
      x=x+d
      write(10,'(3F15.8)')x
    end do
    close(10)
end do
End Program
eof
ifort produce_POSCAR.f90
./a.out
cd ..


#!/bin/bash
rm phvel.dat
q_num=72    #!!!!!!!  No. of irred. q-points
mod=12      #!!!!!!!  No. of phonon branches
for ((i=0;i<${q_num};i=i+1))
do
grep -A $(echo "${mod}+1"|bc) "# Irreducible k point" K3Sb_dos.phvel_all | head -$(echo "$i*($mod+3)+${mod}+2"|bc) | tail -${mod} | awk -v val=$j '{print val, $4,$5}' >> phvel.dat
#!!!!!!!!!!!# the *.phvel_all file name above shoud be changed with system
done

cat > phvel_1.f90 << eof
Program phvel_plot
implicit none
integer:: q_num=${q_num},mod=${mod},i,j
real,allocatable ::frequency(:,:),phvel(:,:)
open(unit=1,File='phvel.dat')
open(unit=3,File='phvel_plot.dat')
allocate(frequency(q_num,mod))
allocate(phvel(q_num,mod))

do i=1,q_num
  do j=1,mod
    read(1,*)frequency(i,j),phvel(i,j)
  end do
end do

do j=1,mod
  do i=1,q_num
    write(3,'(2F20.10)')frequency(i,j),phvel(i,j)
  end do
  write(3,*)
end do

end program
eof
ifort phvel_1.f90
./a.out

cat > phvel_2.f90 << eof
Program phvel_plot
implicit none
integer:: q_num=${q_num},mod=${mod},i,j,k
real,allocatable ::fre_vel(:,:,:)
open(unit=1,File='phvel.dat')
open(unit=3,File='phvel_plot_2.dat')
allocate(fre_vel(q_num,mod,2))

do i=1,q_num
  do j=1,mod
    read(1,*)(fre_vel(i,j,k),k=1,2)
  end do
end do

do i=1,q_num
  write(3,'(50F20.10)')(fre_vel(i,j,:),j=1,mod)
end do

end program
eof
ifort phvel_2.f90
./a.out

rm phvel_1.f90 phvel_2.f90 a.out


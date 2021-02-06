! the ph_linewidth is in unit [cm^(-1)]
! 1/tau=2*(anh_ph_linewidth + iso_selfenergy)
Program three_ph_scattering
implicit none
character*9 line
integer:: q_num,mode_num,i,j,k,dege,pass_1,pass_2
real,allocatable ::fre_scat(:,:,:),iso_fre_scat(:,:,:)
open(unit=1,File='K3BrO_rta.result') ! *******.result file is changed with system
open(unit=2,File='third_ph_sca.dat') ! anh_ph_linewidth
open(unit=3,File='third_ph_sca_2.dat') 
open(unit=4,File='K3BrO_rta.self_isotope') ! *****. changed with system
open(unit=6,File='isot_sca.dat')     ! iso_ph_selfenergy
open(unit=7,File='isot_sca_2.dat')
read(1,*)
read(1,*)
read(1,*)mode_num
mode_num=3*mode_num
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)q_num

allocate(fre_scat(q_num,mode_num,2))
allocate(iso_fre_scat(q_num,mode_num,2))

do while(line.ne.'#K-point ')
  read(1,100)line
  100 format(A9)
  write(*,*)line
end do

do i=1,q_num
  do j=1,mode_num
    read(1,*)pass_1,pass_2,fre_scat(i,j,1)
  end do
end do

read(1,*)
read(1,*)
read(1,*)

do i=1,q_num
  do j=1,mode_num
    read(1,*)
    read(1,*)
    read(1,*)dege
    do k=1,dege
      read(1,*)
    end do
    read(1,*)fre_scat(i,j,2)
    read(1,*)
  end do
end do

do j=1,mode_num
  do i=1,q_num
    write(2,'(2F20.10)')fre_scat(i,j,:)
  end do
  write(2,*)
end do

do i=1,q_num
  write(3,'(50F20.10)')(fre_scat(i,j,:),j=1,mode_num)
end do

read(4,*)
read(4,*)
read(4,*)

do i=1,q_num
 read(4,*)
 read(4,*)
 do j=1,mode_num
  read(4,*)pass_1,pass_2,iso_fre_scat(i,j,:)
 end do
 read(4,*)
end do

do j=1,mode_num
  do i=1,q_num
    write(6,'(2F20.10)')iso_fre_scat(i,j,:)
  end do
  write(6,*)
end do

do i=1,q_num
  write(7,'(50F20.10)')(iso_fre_scat(i,j,:),j=1,mode_num)
end do

end program


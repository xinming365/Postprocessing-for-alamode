Program kappa_spectra
implicit none
integer::i,j,k,row=200
real::pass_t
real,allocatable ::spec(:,:),culm_spec(:,:)
open(unit=1,File='pc_rta.kl_spec')    !!!!!!  change this file for different file
open(unit=2,File='kappa_spec_culm.dat')

write(*,*)'Input row of the *.kl_spec file:'
read(*,*) row
allocate(spec(row,4))
allocate(culm_spec(row,3))

read(1,*)
read(1,*)
do i=1,row
 read(1,*)pass_t,spec(i,:)
end do

culm_spec=0
do j=1,3
 do i=1,row
  do k=1,i
   culm_spec(i,j)=culm_spec(i,j)+spec(k,j+1)
  end do
 end do
end do

culm_spec=culm_spec*(spec(2,1)-spec(1,1))
do i=1,row
 write(2,300,advance='no')spec(i,:),culm_spec(i,:)
 write(2,*)
end do
300 format(7f12.7)
end program


program read_CCOORD
implicit none
integer :: i,ii,iii,j,jj,jjj,defectcnt,noneed,numsteps,numatmtypes,step
integer, allocatable :: countypes(:),indexnum(:)
character(len=2), allocatable :: atomtypes(:)
character(len=2) ::inputypes(1:30)
character (len=30) :: skip,filename
real, allocatable :: fracx(:),fracy(:),fracz(:)
real :: v(1:3,1:3),x,y,z,a,b,c,vol,al,be,ga,time,mass(1:30),amu
amu=0
open(20,file='indvidualcoords.dat')
open(30,file='input')
write(20,'(A20)',advance='no')'Time'
read(30,*)numatmtypes
 do i=1,numatmtypes
  read(30,'(A2,F6.4)')inputypes(i),mass(i)
  write(20,'(8x,A2)',advance='no')inputypes(i)
 enddo
 write(20,*)
 close(30)
open(10,file='CCOORD', status='old')
read(10,*)
read(10,'(A20,I10)')skip,numsteps
do iii=1,numsteps
 read(10,'(A30,I10,I10,f20.6)')skip,defectcnt,step,time
 write(20,'(f20.6)',advance='no')time
 allocate(atomtypes(1:defectcnt))
 allocate(indexnum(1:defectcnt))
 allocate(countypes(1:numatmtypes))
 countypes(:)=0
 do j = 1,3
  read(10,*)(v(j,jj),jj=1,3)
 enddo
 
  do jj=1,defectcnt
   read(10,'(a2,I10)')atomtypes(jj),indexnum(jj)
   do i=1,numatmtypes
    if(inputypes(i).eq.atomtypes(jj))then
     countypes(i)=countypes(i)+1
     
    endif
   enddo
  enddo 
  do i= 1,numatmtypes
   write(20,'(I10)',advance='no')countypes(i)
  enddo
  write(20,*)

  deallocate(atomtypes)
  deallocate(indexnum)
  deallocate(countypes)


enddo
 end

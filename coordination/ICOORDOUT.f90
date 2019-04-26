program read_CCOORD
implicit none
integer :: i,ii,iii,j,jj,jjj,defectcnt,noneed,numsteps,numatmtypes,numatms,steps,&
        numcoorded,inumcoorded,coordination,icoordination,numuniquepairs
character(len=2), allocatable :: atomtypes(:)
character(len=2) ::inputypes(1:30),atm1,atm2,atm1check,atm2check
character (len=40) :: skip,filename
real, allocatable :: fracx(:),fracy(:),fracz(:)
real :: v(1:3,1:3),x,y,z,a,b,c,vol,al,be,ga,testz,testy,testx,mass(1:30),amu,time
Logical :: newatom,itsopen
numsteps=6
amu=0
jj=0
open(30,file='input')
read(30,*)numatmtypes,numatms,numuniquepairs
 do i=1,numatmtypes
 read(30,*)inputypes(i)
 enddo
 close(30)
open(10,file='ICOORD', status='old')
read(10,*)
read(10,*)
do i = 1 ,numatms
read(10,*)
enddo
do ii= 1,numsteps
read(10,'(a40,I10,F20.6)')skip,steps,time
do i = 1,2*numuniquepairs
jj=0
if(i.eq.1)then
read(10,*)atm1,atm2,icoordination,inumcoorded
endif
read(10,*)atm1check,atm2check,coordination,numcoorded
write (filename,'(A1,A1,A1,A4)')atm1,'-',atm2,'.dat'
if(ii.eq.1)then
open (unit=20+i, access='sequential', file=filename)
write(20+i,'(A9,15I7)')'time',0,1,2,3,4,5,6,7,8,9,10,11,12,13,14
end if

 write(20+i,'(F9.6,I7)',advance='no')time,inumcoorded
do while(atm1.eq.atm1check .and. atm2.eq.atm2check)
write(20+i,fmt='(I7)',advance='no')numcoorded
read(10,fmt='(A3)',advance='no')atm1check !,atm2check,coordination,numcoorded




If(atm1.eq.'O' .and. atm2.eq.'Si')then
        print*,coordination,numcoorded,inumcoorded
endif

if(atm1check.eq.'oo')then
        exit
else
read(10,*)atm2check,coordination,numcoorded
endif
enddo
atm1=atm1check
atm2=atm2check
inumcoorded=numcoorded
write(20+i,*)
enddo
enddo



 end
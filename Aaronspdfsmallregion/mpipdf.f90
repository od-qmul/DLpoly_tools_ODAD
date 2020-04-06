Program hellompi
!use mpi
implicit none
integer::ierror,myRank,uniSize,version,subversion,i,ii,iii,iiii,j,jj,jjj,numatoms,numregion,numcom,&
         iMyName,counter,x0,x1,y0,y1,z0,z1,nskip,&
         numinx,numinz,numiny,indx,numinbox,k,kk,kkk,bins(1:200),numbins,&
         biginx,biginy,biginz
real :: delr,bin,rcut,v(1:3,1:3),rarea,&
        numdenbox,pi,volume,rcom,xcom,xbig,ycom,ybig,zcom,zbig,&
        b1,b2,b3,constant
real , allocatable ::x(:),y(:),z(:)
integer , allocatable::boxlist(:,:,:,:),bboxlist(:,:,:,:)
character( len=2) ::atomtype
!character(len=MPI_MAX_PROCESSOR_NAME)::myName
character(len=60)::order

character(len=20) :: x0char,x1char,y0char,y1char,z0char,z1char

x0 = 0.0
x1 = 0.0
y0 = 0.0
y1 = 0.0
z0 = 0.0
z1 = 0.0

if (command_argument_count().NE.6) then
  write(*,*)"Error, 6 command line arguments required, x0, x1, y0, y1, z0, z1"
  stop
end if

call get_command_argument(1,x0char)
call get_command_argument(2,x1char)
call get_command_argument(3,y0char)
call get_command_argument(4,y1char)
call get_command_argument(5,z0char)
call get_command_argument(6,z1char)

read(x0char,*)x0
read(x1char,*)x1
read(y0char,*)y0
read(y1char,*)y1
read(z0char,*)z0
read(z1char,*)z1

write(*,'(6F12.7)') x0,x1,y0,y1,z0,z1


counter=0

!call MPI_Init(ierror)
!call MPI_Comm_size(MPI_COMM_WORLD,uniSize,ierror)
!call MPI_Comm_rank(MPI_COMM_WORLD,myRank, ierror)
!call MPI_Get_processor_name(myName,iMyName, ierror)
!call MPI_Get_version(version,subversion, ierror)
! definig region to search
! 29.246, 33.732, 83.796 
! 10,10,60
!x0=-150.0 ; x1=-100.0; y0=-150.0 ; y1=-100.0; z0=-150.0 ; z1=-100.0

rcut=10.0
bin=0
delr=0.05
bins(1:200)=0
numinx=abs(x1-x0)/rcut;numiny=abs(y1-y0)/rcut;numinz=abs(z1-z0)/rcut
biginx=(abs(x1-x0)-30)/rcut+1 ; biginy=(abs(y1-y0)-30)/rcut+1 ;biginz=(abs(z1-z0)-30)/rcut +1
print*,biginx,biginy,biginz
print*,numinx
numregion=0
numcom=0
numbins=(rcut/delr)
allocate(boxlist(1:numinx,1:numiny,1:numinz,0:10000))
allocate(bboxlist(1:biginx,1:biginy,1:biginz,0:1000000))
boxlist(:,:,:,:)=0
bboxlist(:,:,:,:)=0
if(myRank==0)then
 open(20,file='test')
 open(30,file='pdf.dat')
 open(10,file='REVCON', status='old')
 read(10,*)
 read(10,*)nskip,nskip,numatoms
 
 do j = 1,3
  read(10,*)(v(j,jj),jj=1,3)
 enddo
b1=sqrt(v(1,1)**2 + v(1,2)**2 +v(1,3)**2)
b2=sqrt(v(2,1)**2 + v(2,2)**2 +v(2,3)**2)
b3=sqrt(v(3,1)**2 + v(3,2)**2 +v(3,3)**2) 
volume=b1*b2*b3
numdenbox=numatoms/volume
print*,numdenbox
 allocate(x(1:numatoms))
 allocate(y(1:numatoms))
 allocate(z(1:numatoms))
 counter=0
 do i= 1, numatoms
  read(10,'(A2,I16)')atomtype,indx
  read(10,*)x(i),y(i),z(i)
  read(10,*)
  read(10,*)
 
  do j= 1,numinx
  do jj=1,numiny
  do jjj=1,numinz
   
  if(x(i).ge.(x0+(j-1)*rcut) .and. x(i).lt.(x0+(j)*rcut))then
  if(y(i).ge.(y0+(jj-1)*rcut) .and. y(i).lt.(y0+(jj)*rcut))then
  if(z(i).ge.(z0+(jjj-1)*rcut) .and. z(i).lt.(z0+(jjj)*rcut))then
   numregion=numregion+1
     boxlist(j,jj,jjj,0)=boxlist(j,jj,jjj,0)+1
     
     boxlist(j,jj,jjj,boxlist(j,jj,jjj,0))=indx     
    endif
   endif
  endif
  enddo
  enddo
  enddo
  do j= 1,biginx
  do jj=1,biginy
  do jjj=1,biginz
!   print*,x0+(j-1)*rcut,x0+30+(j-1)*rcut
  if(x(i).ge.(x0+(j-1)*rcut) .and. x(i).lt.(x0+30+(j-1)*rcut))then
  if(y(i).ge.(y0+(jj-1)*rcut) .and. y(i).lt.(y0+30+(jj-1)*rcut))then
  if(z(i).ge.(z0+(jjj-1)*rcut) .and. z(i).lt.(z0+30+(jjj-1)*rcut))then
     
     bboxlist(j,jj,jjj,0)=bboxlist(j,jj,jjj,0)+1
!     print*,bboxlist(j,jj,jjj,0)

     bboxlist(j,jj,jjj,bboxlist(j,jj,jjj,0))=indx
    endif
   endif
  endif

  enddo
  enddo
  enddo
 enddo
 print*,boxlist(2,2,2,0),numregion
bin=0
do i = 10,numbins
  bin=delr*(i-1)
 do jj = 2,numinx-1
 do j = 2,numiny-1  
 do ii= 2,numinz-1
 
    do k = 1,int(boxlist(jj,j,ii,0))
     
     do kk = 1,int(bboxlist(jj-1,j-1,ii-1,0))
      if(boxlist(jj,j,ii,k).ne.bboxlist(jj-1,j-1,ii-1,kk))then     
       xcom=x(boxlist(jj,j,ii,k))      
       ycom=y(boxlist(jj,j,ii,k))
       zcom=z(boxlist(jj,j,ii,k))
       xbig=x(bboxlist(jj-1,j-1,ii-1,kk))
       ybig=y(bboxlist(jj-1,j-1,ii-1,kk))
       zbig=z(bboxlist(jj-1,j-1,ii-1,kk))

 rcom=sqrt( (xcom-xbig)**2 + (ycom-ybig)**2 + (zcom-zbig)**2)
        
         
        if((rcom.ge.(bin)).and.(rcom.lt.(bin+delr)))then
         bins(i)=bins(i)+1
         
!          print*,rcom,xcom,ycom,zcom,xbig,zbig,boxlist(jj,j,ii,k),bboxlist(jj-1,j-1,ii-1,kk)
        endif
        else
       endif
     enddo
    enddo 
   enddo
  enddo
 enddo
enddo

constant=16*atan(1.d0)/3
print*,constant

 do jj = 2,numinx-1
 do j = 2,numiny-1
 do ii= 2,numinz-1

  
      numcom=numcom+int(boxlist(jj,j,ii,0))
enddo
enddo
enddo
print*,numcom
write(30,'("INput coordinates: "(6F12.7))') x0,x1,y0,y1,z0,z1
do i = 1,numbins
 write(30,*)((delr*(i-1)+delr*i))/2,real(bins(i))/(((((delr*i)**3) - (delr*(i-1))**3))*numdenbox*constant*numcom)

!/(16*atan(1.d0)*((delr*(i-1))**2)*numdenbox*counter*0.05))

 enddo



  do j= 1,numinx
  do jj=1,numiny
  do jjj=1,numinz
!   write(20,*)boxlist(j,jj,jjj,0)
   do k = 1,boxlist(j,jj,jjj,0)
!   write(20,*)boxlist(j,jj,jjj,k)
   enddo
  enddo
  enddo
  enddo





 
else
! call MPI_SEND(order,3,MPI_CHARACTER,0,10,MPI_COMM_WORLD,status,ierror)
endif






























end

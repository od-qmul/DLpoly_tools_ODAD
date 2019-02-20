program density_revcon
implicit none 
integer :: densatoms,i,ii,iii,j,jj,jjj,k,kk,kkk,numatoms,ignore,numinx,numiny,numinz,shift
character( len= 8) :: atomtypes(1:3000),atmtypes(1:3000)
character( len= 80) :: line
real, allocatable :: mass(:),x(:),y(:),z(:)
real :: x0,x1,y0,y1,z0,z1,xmid,ymid,zmid,totalmass,v(1:3,1:3),boxlen,density,volume,amu,a,b,c
amu=1.66054e-27
open(10,file='densinput', status='old')
read(10,*)densatoms
!allocate(atomtypes(1:densatoms))
allocate(mass(1:densatoms))
do i = 1, densatoms
read(10,*)atomtypes(i),mass(i)
enddo
close(10)

boxlen=15.0
shift=1
open(20,file='REVCON',status='old')
open(30,file='densmap.dat')
read(20,*)
read(20,*)ignore,ignore,numatoms
!allocate(atmtypes(1:numatoms))
allocate(x(1:numatoms))
allocate(y(1:numatoms))
allocate(z(1:numatoms))
do j = 1,3
read(20,*)(v(j,jj),jj=1,3)

enddo
do i = 1, numatoms
read(20,'(A8)')atmtypes(i)
read(20,*)x(i),y(i),z(i)
read(20,*)
read(20,*)
enddo
a=sqrt(v(1,1)**2 + v(1,2)**2 + v(1,3)**2)
b=sqrt(v(2,1)**2 + v(2,2)**2 + v(2,3)**2)
c=sqrt(v(3,1)**2 + v(3,2)**2 + v(3,3)**2)
numinx=(a-boxlen)/shift
numiny=(b-boxlen)/shift
numinz=(c-boxlen)/shift
write(30,*)'number of slices =',numinz,numinx,numiny
 do kk=1,numinz
   write(30,*)'slice',kk
 do jj=1,numiny
    
  do ii=1,numinx 
    totalmass=0
     x0=(-(a/2.0)+ii*shift)
     y0=(-(b/2.0)+jj*shift)
     z0=(-(c/2.0)+kk*shift)
!     print*,z0,z1,-c/2.0
     x1=(-(a/2.0)+boxlen+ii*shift)
     y1=(-(b/2.0)+boxlen+jj*shift)
     z1=(-(c/2.0)+boxlen+kk*shift)
     xmid=X0+boxlen/2
     ymid=y0+boxlen/2
     zmid=z0+boxlen/2
     do i=1,numatoms
     if(x(i).ge.x0 .and. x(i).le.x1)then
      if(y(i).ge.y0 .and. y(i).le.y1)then
       if(z(i).ge.z0 .and. z(i).le.z1)then
        do jjj = 1, densatoms
         if(atmtypes(i).eq.atomtypes(jjj))then         
          totalmass=totalmass+mass(jjj)
         endif
        enddo 
       endif
      endif  
     endif 
     enddo
     volume=(boxlen)**3
     
     density=(totalmass*amu/volume)*(1E3/(1E-8)**3)
     write(30,'(F7.5,1x)',advance='no')density
  enddo
     write(30,*)
 enddo
  
 enddo
end

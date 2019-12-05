program find_defects_from_ICOORD
integer :: numatoms,numpairs,perfectcoord(1:100),indx,coord,coordpair,&
           numdefects,i,ii,iii,j,jj,jjj,Next,aux,tempindx(1:100),&
           numinx,numiny,numinz,shift,numregion,numuniq,numtest,&
           ndefs,numbins,bins(1:200),bin,numdefs,testdefs
real    :: v(1:3,1:3),a,b,c,amu,denconstant,delr,r
real ,allocatable :: atomass(:),boxmass(:,:,:),x1(:),&
                     y1(:),z1(:),x2(:),y2(:),z2(:)
integer ,allocatable :: defectindx1(:),defectindx2(:),truedefects(:),defects(:)
character(len=2)::coordpairs(1:100,1:2),atm1,atm2,currentatm,tempatm(1:100)
character(len=2), allocatable::atm(:),atomtypes(:)
character(len=80):: title
logical:: correctatm,found
boxlen=20
shift=20
amu=1.66054e-27
numtest=0
numdefects=0

delr=0.001
bins=0
open(100,file='test')
open(10,file='icoordinput',status='old')
read(10,*)numpairs,numatoms,numuniq
do i = 1,numpairs
read(10,*)coordpairs(i,1),coordpairs(i,2),perfectcoord(i)
enddo
allocate(atomtypes(1:numuniq))
allocate(atomass(1:numuniq))
do i=1,numuniq
read(10,*)atomtypes(i),atomass(i)
enddo
allocate(truedefects(1:numatoms))
allocate(defectindx1(1:numatoms))
allocate(defectindx2(1:numatoms))
allocate(atm(1:numatoms))
allocate(x1(1:numatoms))
allocate(y1(1:numatoms))
allocate(z1(1:numatoms))
allocate(x2(1:numatoms))
allocate(y2(1:numatoms))
allocate(z2(1:numatoms))

! reading ICOORD indexes to identify defects
currentatm=''
numdefects=0
numberdefects2=0

open(30,file='CONFIG1',status='old')
read(30,'(a80)')title
read(30,*)
do j = 1,3
  read(30,*)(v(j,jj),jj=1,3)
 enddo
a=sqrt(v(1,1)**2 + v(1,2)**2 + v(1,3)**2)
b=sqrt(v(2,1)**2 + v(2,2)**2 + v(2,3)**2)
c=sqrt(v(3,1)**2 + v(3,2)**2 + v(3,3)**2)
numinx=nint((a-boxlen)/shift +1)
numiny=nint((b-boxlen)/shift +1)
numinz=nint((c-boxlen)/shift +1)


do i = 1,numatoms

read(30,*)atm(i),indx
read(30,*)x1(i),y1(i),z1(i)
read(30,*)
read(30,*)
enddo
close(30)
print*,'read config'





open(20, file='ICOORD1',status='old')
read(20,*)
read(20,*)
do i = 1, numatoms
   read(20,Fmt='(i12,1x,i12)')indx,coord
   write(100,*)indx,coord,i
 
   if(coord.gt.0)then 
    if(currentatm.ne.atm(indx))then
     do ii= 1,numpairs
        if(coordpairs(ii,1).eq.atm(indx))then
          currentatm=atm1
          coordpair=ii
           if(perfectcoord(coordpair).ne.coord)then
             numdefects=numdefects+1
             defectindx1(numdefects)=indx
           else
           endif
        else  
       endif
     enddo

    else
        if(perfectcoord(coordpair).ne.coord)then
         numdefects=numdefects+1
         defectindx1(numdefects)=indx 
        else
        endif
    endif
   else

   endif 
enddo
print*,numdefects 


close(20)
close(10)




open(40,file='CONFIG2',status='old')
read(40,'(a80)')title
read(40,*)
do j = 1,3
  read(40,*)(v(j,jj),jj=1,3)
 enddo
a=sqrt(v(1,1)**2 + v(1,2)**2 + v(1,3)**2)
b=sqrt(v(2,1)**2 + v(2,2)**2 + v(2,3)**2)
c=sqrt(v(3,1)**2 + v(3,2)**2 + v(3,3)**2)
numinx=nint((a-boxlen)/shift +1)
numiny=nint((b-boxlen)/shift +1)
numinz=nint((c-boxlen)/shift +1)


do i = 1,numatoms

read(40,*)atm(i),indx
read(40,*)x2(i),y2(i),z2(i)
read(40,*)
read(40,*)
enddo
close(40)
print*,'read config'


open(50, file='ICOORD2',status='old')
read(50,*)
read(50,*)
do i = 1, numatoms
   read(50,Fmt='(i12,1x,A8,i12)')indx,atm1,coord
!   write(100,*)indx,atm1,coord
   if(coord.gt.0)then
    if(currentatm.ne.atm(indx))then
     do ii= 1,numpairs
        if(coordpairs(ii,1).eq.atm(indx))then
          currentatm=atm1
          coordpair=ii
           if(perfectcoord(coordpair).ne.coord)then
             numdefects2=numdefects2+1
             defectindx2(numdefects2)=indx
           else
             
           endif
        else
     
        endif
     enddo

    else
        if(perfectcoord(coordpair).ne.coord)then
         numdefects2=numdefects2+1
         defectindx2(numdefects2)=indx
            
        else
        
        endif
    endif
   else
   
   endif
enddo
print*,numdefects2
numdefs=0


do ii = 1,numdefects2
found=.false.
 do i =1,numdefects
    if(defectindx2(ii).eq.defectindx1(i))then
      found=.true.
    endif
  enddo
    if(found.eqv. .false.)then
      numdefs=numdefs+1
      truedefects(numdefs)=defectindx2(ii)
    endif

enddo
print*,numdefs
allocate(defects(1:numdefs))
ndefs=0
do i = 1,numdefs
x2(truedefects(i))=x2(truedefects(i)) -V(1,1)*anint(x2(truedefects(i))/v(1,1))
y2(truedefects(i))=y2(truedefects(i)) -V(2,2)*anint(y2(truedefects(i))/v(2,2))
z1(truedefects(i))=z2(truedefects(i)) -V(3,3)*anint(z2(truedefects(i))/v(3,3))
x1(truedefects(i))=x1(truedefects(i)) -V(1,1)*anint(x1(truedefects(i))/v(1,1))
y1(truedefects(i))=y1(truedefects(i)) -V(2,2)*anint(y1(truedefects(i))/v(2,2))
z1(truedefects(i))=z1(truedefects(i)) -V(3,3)*anint(z1(truedefects(i))/v(3,3))




r=( ( x2(truedefects(i)) - x1(truedefects(i)))**2 + &
( y2(truedefects(i)) - y1(truedefects(i)))**2 + &
( z2(truedefects(i)) - z1(truedefects(i)))**2 )

if(r.gt.0.75)then
ndefs=ndefs+1
defects(ndefs)=truedefects(i)
endif
enddo




open(110,file='defects.cfg')
    write(110,fmt='(A21,1x,I7)') 'Number of particles =',ndefs
    write(110,fmt='(A3,1x,F5.1,1x,A8)') 'A =',1.0,'Angstrom'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(1,1) =',v(1,1),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(1,2) =',v(1,2),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(1,3) =',v(1,3),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(2,1) =',v(2,1),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(2,2) =',v(2,2),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(2,3) =',v(2,3),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(3,1) =',v(3,1),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(3,2) =',v(3,2),'A'
    write(110,fmt='(A9,1x,F7.2,1x,A1)') 'H0(3,3) =',v(3,3),'A'


do i =1,ndefs
write(110,fmt='(F6.2,1x,A2,6F9.4)')atomass(1),atm(truedefects(i)),&
x2(defects(i))/v(1,1),y2(defects(i))/v(2,2),z2(defects(i))/v(3,3),&
0.d0,0.d0,0.d0

enddo

end

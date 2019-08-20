program find_defects_from_ICOORD
integer :: numatoms,numpairs,perfectcoord(1:100),indx,coord,coordpair,&
           numdefects,i,ii,iii,j,jj,jjj,Next,aux,tempindx(1:100),&
           numinx,numiny,numinz,shift,numregion,numuniq,numtest,&
           ndef,numbins,bins(1:200),bin
real    :: v(1:3,1:3),a,b,c,amu,denconstant,delr
real ,allocatable :: defects(:,:),atomass(:),boxmass(:,:,:),x(:),y(:),z(:)
integer ,allocatable :: defectindx(:),boxlist(:,:,:,:)
character(len=2)::coordpairs(1:100,1:2),atm1,atm2,currentatm,tempatm(1:100)
character(len=2), allocatable::atm(:),atomtypes(:)
character(len=80):: title
logical:: correctatm
boxlen=10
shift=10
amu=1.66054e-27
numtest=0
numdefects=0

delr=0.001
bins=0
open(120,file='denshisto.dat')
open(100,file='densmap.dat')
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

allocate(defectindx(1:numatoms))
allocate(defects(1:numatoms,1:3))
allocate(atm(1:numatoms))
allocate(x(1:numatoms))
allocate(y(1:numatoms))
allocate(z(1:numatoms))

! reading ICOORD indexes to identify defects
currentatm=''
numdefects=0

open(30,file='CONFIG',status='old')
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
read(30,*)x(i),y(i),z(i)
read(30,*)
read(30,*)
enddo
close(30)
print*,'read config'





open(20, file='ICOORD',status='old')
read(20,*)
read(20,*)
do i = 1, numatoms
   read(20,Fmt='(i12,1x,i12)',advance='no')indx,coord
!   write(100,*)indx,coord
   if(coord.gt.0)then 
    if(currentatm.ne.atm(indx))then
     do ii= 1,numpairs
        if(coordpairs(ii,1).eq.atm(indx))then
          currentatm=atm1
          coordpair=ii
           if(perfectcoord(coordpair).ne.coord)then
             numdefects=numdefects+1
             defectindx(numdefects)=indx
             read(20,*)(tempindx(iii),iii=1,coord),(tempatm(jj),jj=1,coord)
              do jjj= 1,coord
                if(tempatm(jjj).eq.coordpairs(ii,2))then
                   numdefects=numdefects+1
                   defectindx(numdefects)=tempindx(jjj)
                endif
              enddo
           else
             read(20,*)
           endif
        else
         read(20,*)       
        endif
     enddo

    else
        if(perfectcoord(coordpair).ne.coord)then
         numdefects=numdefects+1
         defectindx(numdefects)=indx 
             read(20,*)(tempindx(iii),iii=1,coord),(tempatm(jj),jj=1,coord)
                        
          do jjj = 1, coord
             
            if(tempatm(jjj).eq.coordpairs(coordpair,2))then
               
               numdefects=numdefects+1
               defectindx(numdefects)=tempindx(jjj)
             
            endif 
          enddo
        else
         read(20,*)
        endif
    endif
   else
   read(20,*)
   endif 
enddo
print*,numdefects 
!do i= 1,numdefects
! write(100,*)
!enddo





close(20)
close(10)
!sorting defectindx



!open(40,file='icoordefects.xyz')
!write(40,*)numdefects
!write(40,*)title
x0=-a/2
y0=-b/2
z0=-c/2
allocate(boxlist(1:numinx,1:numiny,1:numinz,0:100000))
allocate(boxmass(1:numinx,1:numiny,1:numinz))
boxlist(:,:,:,:)=0
boxmass(:,:,:)=0
print*,'Number of slices=',numinx,numiny,numinz

do i= 1,numdefects
 do j= 1,numinx
  do jj=1,numiny
   do jjj=1,numinz
    if(x(defectindx(i)).ge.(x0+(j-1)*shift) .and. x(defectindx(i)).lt.(x0+(j)*shift))then
     if(y(defectindx(i)).ge.(y0+(jj-1)*shift) .and. y(defectindx(i)).lt.(y0+(jj)*shift))then
      if(z(defectindx(i)).ge.(z0+(jjj-1)*shift) .and. z(defectindx(i)).lt.(z0+(jjj)*shift))then
        boxlist(j,jj,jjj,0)=boxlist(j,jj,jjj,0)+1
        boxlist(j,jj,jjj,boxlist(j,jj,jjj,0))=defectindx(i)
        numtest=numtest+1
!        write(100,*)numtest
        do ii=1,numuniq
        if(atm(defectindx(i)).eq.atomtypes(ii))then
        boxmass(j,jj,jjj)=boxmass(j,jj,jjj)+atomass(ii)
        endif
        enddo
      endif
     endif
    endif
   enddo
  enddo
 enddo
enddo

denconstant=(amu/(boxlen)**3)*(1E3/(1E-8)**3)
write(100,*)numinx
write(100,*)numiny
do jjj= 1,numinz
write(100,*)'Slice',(z0+(jjj-1)*shift),'to',(z0+(jjj)*shift)
do jj=1,numiny
do j=1,numinx
write(100,'(F7.5,1x)',advance='no')boxlist(j,jj,jjj,0)/(boxlen)**3
!boxmass(j,jj,jjj)*denconstant
enddo
write(100,*)
enddo
enddo
open(110,file='defects.xyz')
write(110,*)numdefects
write(110,*)title
do i =1,numdefects
write(110,*)atm(defectindx(i)),x(defectindx(i)),y(defectindx(i)),z(defectindx(i))
enddo
bin=0
numbins=50
do i = 1,numbins
 bin=(i-1)
! print*,bin
 do jjj=1,numinz
  do jj=1,numiny
   do j=1,numinx
       if(bin.eq.boxlist(j,jj,jjj,0))then
!          write(120,*)bin,boxlist(j,jj,jjj,0)
          bins(i)=bins(i)+1
       endif
   enddo
  enddo
 enddo
 write(120,*)bin,bins(i)
enddo
















end

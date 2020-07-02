program find_defects_from_ICOORD
integer :: numatoms,numpairs,perfectcoord(1:100),coord,coordpair,&
           numdefects,i,ii,iii,j,jj,jjj,Next,aux,tempindx(1:100),&
           numinx,numiny,numinz,shift,numregion,numuniq,numtest,&
           ndef,numbins,bin
real    :: v(1:3,1:3),a,b,c,amu,denconstant,delr,coordaverage,binmin
real ,allocatable :: defects(:,:),atomass(:),boxmass(:,:,:),&
                    x(:),y(:),z(:),bins(:)
integer ,allocatable :: indx(:),boxlist(:,:,:,:),coordnum(:)
character(len=2)::coordpairs(1:100,1:2),atm1,atm2,currentatm,tempatm(1:100)
character(len=2), allocatable::atm(:),atomtypes(:)
character(len=80):: title
logical:: correctatm,file_exists
boxlen=10
shift=10
amu=1.66054e-27
numtest=0
numdefects=0

delr=0.01
numbins=1.0/0.01

allocate(bins(1:numbins))
bins=0
open(50,file='test')
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

allocate(indx(1:numatoms))
allocate(coordnum(1:numatoms))
allocate(atm(1:numatoms))
allocate(x(1:numatoms))
allocate(y(1:numatoms))
allocate(z(1:numatoms))

! reading ICOORD indexes to identify defects
currentatm=''
numdefects=0

inquire(file="REVCON", exist=file_exists)

if (file_exists .eqv. .true.) then
  open(30,file='REVCON',status='old')
else
  open(30,file='CONFIG',status='old')
end if

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


x0=-a/2
y0=-b/2
z0=-c/2
allocate(boxlist(1:numinx,1:numiny,1:numinz,0:10000))
allocate(boxmass(1:numinx,1:numiny,1:numinz))
boxlist(:,:,:,:)=0
boxmass(:,:,:)=0
print*,'Number of slices=',numinx,numiny,numinz

do i = 1,numatoms

read(30,*)atm(i),indx(i)
read(30,*)x(i),y(i),z(i)
read(30,*)
read(30,*)

do j= 1,numinx
do jj=1,numiny
do jjj=1,numinz
if(x(i).ge.(x0+(j-1)*shift) .and. x(i).lt.(x0+(j)*shift))then
if(y(i).ge.(y0+(jj-1)*shift) .and. y(i).lt.(y0+(jj)*shift))then
if(z(i).ge.(z0+(jjj-1)*shift) .and. z(i).lt.(z0+(jjj)*shift))then
boxlist(j,jj,jjj,0)=boxlist(j,jj,jjj,0)+1
boxlist(j,jj,jjj,boxlist(j,jj,jjj,0))=indx(i)
do ii=1,numuniq
if(atomtypes(ii).eq.atm(i))then
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
close(30)
print*,'read config'


close(10)
!sorting defectindx



!open(40,file='icoordefects.xyz')
!write(40,*)numdefects
!write(40,*)title
!x0=-a/2
!y0=-b/2
!z0=-c/2
!allocate(boxlist(1:numinx,1:numiny,1:numinz,0:100000))
!allocate(boxmass(1:numinx,1:numiny,1:numinz))
!boxlist(:,:,:,:)=0
!boxmass(:,:,:)=0
print*,'Number of slices=',numinx,numiny,numinz
     


denconstant=(amu/(boxlen)**3)*(1E3/(1E-8)**3)
write(100,*)numinx-2
write(100,*)numiny-2
do jjj= 2,numinz-1
write(100,*)'Slice',(z0+(jjj-1)*shift),'to',(z0+(jjj)*shift)
do jj=2,numiny-1
do j=2,numinx-1

write(100,'(F7.5,1x)',advance='no') real(boxmass(j,jj,jjj)*denconstant)
!boxmass(j,jj,jjj)*denconstant
enddo
write(100,*)
enddo
enddo



open(110,file='coordhisto.dat')
binmin=6
do i=1,numbins
binmin=6+(i-1)*delr
do jjj=1,numinz
 do jj=1,numiny
  do j=1,numinx
   coordaverage =real(boxlist(j,jj,jjj,1))/real(boxlist(j,jj,jjj,0))
    if(coordaverage  .ge.binmin .and.  coordaverage.le.binmin+delr)then
     bins(i)=bins(i)+1
!      print*,binmin,coordaverage
    endif
   enddo
  enddo
 enddo
 write(110,*)(2*binmin+delr)/2,bins(i)
enddo








end

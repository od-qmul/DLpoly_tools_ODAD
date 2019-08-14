program find_defects_from_ICOORD
integer :: numatoms,numpairs,perfectcoord(1:100),indx,coord,coordpair,&
           numdefects,i,ii,iii,j,jj,jjj,Next,aux,tempindx(1:100),&
           numinx,numiny,numinz,shift,numregion,numuniq,numtest
real    :: v(1:3,1:3),a,b,c,amu,denconstant
real ,allocatable :: defects(:,:),atomass(:),boxmass(:,:,:)
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

open(100,file='testing.dat')
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
! reading ICOORD indexes to identify defects
currentatm=''
numdefects=0
open(20, file='ICOORD',status='old')
read(20,*)
read(20,*)
do i = 1, numatoms
   read(20,Fmt='(i12,1x,a8,i12,1x)',advance='no')indx,atm1,coord
    
    if(currentatm.ne.atm1)then
      
     do ii= 1,numpairs
        if(coordpairs(ii,1).eq.atm1)then
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
enddo
print*,numdefects 
!do i= 1,numdefects
! write(100,*)
!enddo





close(20)
close(10)
!sorting defectindx
print*,'sorting'
do j = numdefects/2, 1, -1
  
    Next =defectindx(j)
    i=j+j
!    print*,Next
   
   do while (i<=numdefects)
     if (i<numdefects)then
       if(defectindx(i)<defectindx(i+1)) i=i+1
     endif
     if(defectindx(i)<=Next)exit
     defectindx(i/2)=defectindx(i) 
     i=i+i
   enddo
  defectindx(i/2)=Next
enddo


   
do j= numdefects, 2, -1
    aux=defectindx(1)
    defectindx(1)=defectindx(j)
    defectindx(j)=aux
    Next=defectindx(1)
    i=2
 
   do while(i<=j-1)
    if(i<j-1)then 
      if(defectindx(i)<defectindx(i+1))i=i+1    
    endif
    if(defectindx(i)<=Next)exit
    defectindx(i/2)=defectindx(i)
    i=i+i
  enddo
  defectindx(i/2)=Next
 enddo
k=1
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
read(30,*)x,y,z
read(30,*)
read(30,*)
   if(indx==defectindx(k))then
   defects(k,1)=x
   defects(k,2)=y
   defects(k,3)=z
   k=k+1
   endif 
enddo
close(30)
open(40,file='icoordefects.xyz')
!write(40,*)numdefects
!write(40,*)title
x0=-a/2
y0=-b/2
z0=-c/2
allocate(boxlist(1:numinx,1:numiny,1:numinz,0:10000))
allocate(boxmass(1:numinx,1:numiny,1:numinz))
boxlist(:,:,:,:)=0
boxmass(:,:,:)=0
print*,'Number of slices=',numinx,numiny,numinz

do i= 1,numdefects
 do j= 1,numinx
  do jj=1,numiny
   do jjj=1,numinz
    if(defects(i,1).ge.(x0+(j-1)*shift) .and. defects(i,1).lt.(x0+(j)*shift))then
     if(defects(i,2).ge.(y0+(jj-1)*shift) .and. defects(i,2).lt.(y0+(jj)*shift))then
      if(defects(i,3).ge.(z0+(jjj-1)*shift) .and. defects(i,3).lt.(z0+(jjj)*shift))then
        boxlist(j,jj,jjj,0)=boxlist(j,jj,jjj,0)+1
        boxlist(j,jj,jjj,boxlist(j,jj,jjj,0))=defectindx(i)
        numtest=numtest+1
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

do jjj= 1,numinz
write(100,*)(z0+(jjj-1)*shift)
do jj=1,numiny
do j=1,numinx
write(100,'(F7.5,1x)',advance='no')boxmass(j,jj,jjj)*denconstant
enddo
enddo
write(100,*)
enddo





end

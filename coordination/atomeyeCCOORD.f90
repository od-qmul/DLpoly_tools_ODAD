program read_CCOORD
implicit none
integer :: i,ii,iii,j,jj,jjj,defectcnt,noneed,numsteps,numatmtypes
character(len=2), allocatable :: atomtypes(:)
character(len=2) ::inputypes(1:30)
character (len=30) :: skip,filename
real, allocatable :: fracx(:),fracy(:),fracz(:)
real :: v(1:3,1:3),x,y,z,a,b,c,vol,al,be,ga,testz,testy,testx,mass(1:30),amu
numsteps=500
amu=0
open(30,file='input')
read(30,*)numatmtypes
 do i=1,numatmtypes
 read(30,*)inputypes(i),mass(i)
 enddo
 close(30)
open(10,file='CCOORD', status='old')
read(10,*)skip
print*,skip
read(10,*)skip
do iii=1,numsteps
read(10,'(A30,I10)')skip,defectcnt
print*,defectcnt
do j = 1,3
  read(10,*)(v(j,jj),jj=1,3)
   enddo

a=sqrt(v(1,1)**2 +v(1,2)**2 +v(1,3)**2)
b=sqrt(v(2,1)**2 +v(2,2)**2 +v(2,3)**2)
c=sqrt(v(3,1)**2 +v(3,2)**2 +v(3,3)**2)
al=acos((v(3,1)*v(2,1) + v(3,2)*v(2,2) + v(3,3)*v(2,3))/abs(c*b))
be=acos((v(3,1)*v(1,1) + v(3,2)*v(1,2) + v(3,3)*v(1,3))/abs(c*a))
ga=acos((v(1,1)*v(2,1) + v(1,2)*v(2,2) + v(1,3)*v(2,3))/abs(a*b))
vol=a*b*c*(sqrt(1-(cos(al)**2) -(cos(be)**2) -(cos(ga)**2) &
        + 2*cos(al)*cos(be)*cos(ga)))
 allocate(atomtypes(1:defectcnt))
 allocate(fracx(1:defectcnt))
 allocate(fracy(1:defectcnt))
 allocate(fracz(1:defectcnt))
  do i =1,defectcnt
   read(10,*)atomtypes(i),noneed,x,y,z
      fracz(i)=(z*a*b*sin(ga))/vol
      fracy(i)=y/b*sin(ga) + (a*c*cos(be)*cos(ga)-cos(al))/vol*sin(ga)
      fracx(i)=x/a - (cos(ga)*y)/a*sin(ga) + (b*c*cos(al)*cos(ga) -cos(be)*z)/vol*sin(ga)
  enddo
 write (filename, '(I0,a4)') iii,'.cfg'
  open (unit=20, access='sequential', file=filename)
          write(20,fmt='(A21,1x,I7)') 'Number of particles =',defectcnt
          write(20,fmt='(A3,1x,F5.1,1x,A8)') 'A =',1.0,'Angstrom'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(1,1) =',v(1,1),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(1,2) =',v(1,2),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(1,3) =',v(1,3),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(2,1) =',v(2,1),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(2,2) =',v(2,2),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(2,3) =',v(2,3),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(3,1) =',v(3,1),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(3,2) =',v(3,2),'A'
          write(20,fmt='(A9,1x,F7.2,1x,A1)') 'H0(3,3) =',v(3,3),'A'
      
      do i =1,defectcnt
       do ii=1,numatmtypes
        If(atomtypes(i).eq.inputypes(ii))then
                     amu=mass(ii)
        endif
       enddo        
      write(20,fmt='(F6.2,1x,A2,6F9.4)')amu,atomtypes(i),fracx(i),fracy(i),fracz(i),&
              0.d0,0.d0,0.d0
      enddo        



deallocate(atomtypes)
deallocate(fracx)
deallocate(fracz)
deallocate(fracy)
enddo



 end

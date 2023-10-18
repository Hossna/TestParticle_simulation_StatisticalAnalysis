      include 'values.mdl'
!=======================================================================
      program momentsDist
      use VALUES
      implicit none

!  this program reads in a set of distribution functions, and then calculates
! various quantities from them.  The distribution functions must be called
! distfunct*.out, where * increases from 000001 in steps of 1.  The x
! position corresponding to the first distribution function
! (distfunct000001.out), and the step by which the x position changes from
! one distribution function to the next are specified in the file
! readDistValues.dat.   
! This program  prints the results in the file momentsdist.out.

!  variables
! nCells: the number of cells contained in the distribution function file
! nVertices: the number of vertices contained in the distribution function
!            file.
! f: used to hold the function needed to be integrated at each of the
!    eight vertices of the cells to be integrated over.
! mom: the moment calculated from the distribution function.
! inputFile: names of the distribution functions being read in.
! numFiles: the number of files containing distribution function.
! lunit: true or false whether fluxes are computed along given unit vertors
!        or not.
! uxyz: coordinates of the unit vector along which fluxes are computed.
      integer i
      integer numFiles
      real*8 uxyz(3),mom
      logical lunit

      namelist /readDistValues/ mass,q,type_integration,numFiles,lunit

!  1.  read in the values from the namelist
      open(unit=9,file='momentsdist.dat',status='old')
      read(9,nml=readDistValues)
      close(9)
!  1.1 open the file containing the unit vectors along which we need to
!      compute the fluxes
      if(lunit) open(unit=9,file='momentsunit.dat',status='old')

!  2.  open a file where the output data is to be printed
      open(unit=15,file='momentsdist.out',status='unknown')

!  3.  do calculations for each of the files, one by one.
      do i=1,numFiles
!  3.1 read in txyzf, and tindex from the file.
         call readFile(i)
         if(lunit) read(9,*)uxyz
   
!  3.2 compute the moment
         call moment(lunit,uxyz,mom)

!  3.3 print the result
         write(15,104)mom,uxyz

!  3.4 deallocate the txyzf, and tindex arrays
         deallocate(txyzf)
         deallocate(tindex)
      enddo

!rm print*,'average=',average
102 format(7es16.6)
103 format(99es14.6)
104 format(es15.5,3f10.4)
      stop
      end

!===========================================================================
      subroutine readFile(fileNum)
      use VALUES
      implicit none

!  this subroutine reads in values for txyzf, and tindex from a file.

!  arguments
!   fileNum: the number of the file to be read in.
      integer fileNum

!  local variables
      integer i,ii,jj
!   line: used when reading in the distribution function
      character(len=132)line
!   inputFile: used to hold the name of the distribution function being
!              read in.
      character(len=60) inputFile


!  1.  make a call to charCompose to get the name of the file to
!      be read in.
      call charCompose('distfunct',fileNum,'.out',inputFile)
      open(unit=2,file=inputFile,status='old')

!  2.  read in the values for txyzf, and tindex from the file.
1     continue
      read(2,101)line
101   format(a)
      if(index(line,'ZONE T=')==0)go to 1
      ii=index(line,'N=')
      jj=index(line,'E=')
      read(line(ii+2:jj),*)nVertices
      read(line(jj+2:132),*)nCells
      allocate(txyzf(4,nVertices))
      allocate(tindex(8,nCells))
      do i=1,nVertices
         read(2,*)txyzf(1,i),txyzf(2,i),txyzf(3,i),txyzf(4,i)
      enddo
      do i=1,nCells
         read(2,*)tindex(1,i),tindex(2,i),tindex(4,i),tindex(3,i),&
              tindex(5,i),tindex(6,i),tindex(8,i),tindex(7,i)
      enddo

      close(2)
!  3.  the values of txyzf were normalized when they were printed
!      in the distribution file, so that the files could be viewed using
!      Vu (when they aren't normalized, the numbers are too small).
!      The values are now converted back to what they should be.
if(1 == 2) then !This was used in the original (Frances's) version
                !of the code, but it is no longer needed.
      if(type_integration==1)then
         do i=1,nVertices
            txyzf(1,i)=txyzf(1,i)*sqrt(2.*mass*q)
            txyzf(2,i)=txyzf(2,i)*sqrt(2.*mass*q)
            txyzf(3,i)=txyzf(3,i)*sqrt(2.*mass*q)
            txyzf(4,i)=txyzf(4,i)/sqrt(2.*mass*q)**3
         enddo
      else
         do i=1,nVertices
            txyzf(1,i)=txyzf(1,i)*sqrt(2.*q/mass)
            txyzf(2,i)=txyzf(2,i)*sqrt(2.*q/mass)
            txyzf(3,i)=txyzf(3,i)*sqrt(2.*q/mass)
            txyzf(4,i)=txyzf(4,i)/sqrt(2.*q/mass)**3
         enddo
      endif
endif
    
      return
      end

!=======================================================================
      subroutine moment(lunit,uxyz,mom)
      use VALUES
      implicit none

!  this subroutine calculates various quantities from the distribution
!  functions.


!  arguments
!   type specifies the quantity to be calculated.
!   calcValue: to contain the quantity calculated.
      real*8 mom,uxyz(3)
      logical lunit

!  local variables
      integer i,j
      real*8 f(8),integral,dotp

!  2.  Compute the flux along unit vector uxyz
!  2.1 if type==1, calculate the particle flux
      if(lunit) then
        do i=1,nCells
          do j=1,8
            dotp=dot_product(txyzf(:,tindex(j,i)),uxyz(:))
            if(dotp <0.d0) then
              f(j)=txyzf(4,tindex(j,i))*abs(dotp)
            else
              f(j)=0.d0
            endif
          enddo
          call integrateCell(i,f,integral)
          mom=mom+integral
        enddo
      endif

      return
      end
!=======================================================================
      subroutine charCompose(ch1,i,ch2,outChar)
      IMPLICIT NONE
!  This returns a character made of ch1, followed by integer i,
!  followed by ch2. This is used to construct the outputFileName.
!
!  arguments
      CHARACTER (LEN=*) :: ch1,ch2,outChar
      INTEGER i
!
!  computation
      if(i < 10) then
        write(outChar,101)ch1//'00000',i,ch2
 101    format(a,i1,a)
      elseif(i < 100) then
        write(outChar,102)ch1//'0000',i,ch2
 102    format(a,i2,a)
      elseif(i < 1000) then
        write(outChar,103)ch1//'000',i,ch2
 103    format(a,i3,a)
      elseif(i < 10000) then
        write(outChar,104)ch1//'00',i,ch2
 104    format(a,i4,a)
      elseif(i < 100000) then
        write(outChar,105)ch1//'0',i,ch2
 105    format(a,i5,a)
      elseif(i < 1000000) then
        write(outChar,106)ch1,i,ch2
 106    format(a,i6,a)
      endif
      return
      end
!=====================================================================
      subroutine integrateCell(index,f,integral)
      use VALUES
      implicit none

!  this subroutine calculates the integral over one cell.

!  arguments
!   index gives the index of the cell over which the integral is
!   to be calculated.
      integer index
!   f holds the function values at the eight vertices
      real*8 f(8)
!   integral is the value of the integral in the cell
      real*8 integral

!  local variables
!   x,y,z: x(1) is the minimum value of the x-momentum in the cell,
!          x(2) is the maximum value of the x-momentum in the cell.
!          Similarly for y, and z.
      real*8 x(2),y(2),z(2)


!  1.  Find the minimum and maximum values of the momentum
!      for the cell.
!   1.1  Minimum values:
      x(1)=txyzf(1,tindex(1,index))
      y(1)=txyzf(2,tindex(1,index))
      z(1)=txyzf(3,tindex(1,index))
!   1.2  Maximum values:
      x(2)=txyzf(1,tindex(2,index))
      y(2)=txyzf(2,tindex(3,index))
      z(2)=txyzf(3,tindex(5,index))

!  2.  calculate the integral.
      integral=0.125*(x(2)-x(1))*(y(2)-y(1))*(z(2)-z(1))*&
           (f(1)+f(2)+f(3)+f(4)+f(5)+f(6)+f(7)+f(8))

!  3.  return
      return
      end

      include 'values.mdl'
!============================================================================
      program readDist
      use VALUES
      implicit none

!  this program reads in a set of distribution functions, and then calculates
!  various quantities from them.  The distribution functions must be called
!  distfunct*.out, where * increases from 000001 in steps of 1.  The x
!  position corresponding to the first distribution function
!  (distfunct000001.out), and the step by which the x position changes from
!  one distribution function to the next are specified in the file
!  readDistValues.dat.   
!  This program  prints the results in the file properties.out, where the 
!  first column is the x positions, the second column is the density, the 
!  third column is the energy, the fourth column is the energy divided by 
!  the density, and the last three columns are the average momentum or 
!  velocity.  

!  variables
      integer i,j,k
!   nCells: the number of cells contained in the distribution function file
!   nVertices: the number of vertices contained in the distribution function
!              file.
!      integer nCells, nVertices
!   f: used to hold the function needed to be integrated at each of the
!      eight vertices of the cells to be integrated over.
      real*8 f(8)
!   integral: used to hold the value of the integral of the function f
!             over a cell.
      real*8 integral
!   density: the density calculated from the distribution function.
!   vx_bar,vy_bar,vz_bar: the average velocity calculated from the distribution
!                         function.
!   px_bar,py_bar,pz_bar: the average momentum calculated from the distribution
!                         function
      real*8 density,vx_bar,vy_bar,vz_bar,px_bar,py_bar,pz_bar,Ti
      real*8 momentum(3),stressTensor(6),pintegr,dnintegr

!   inputFile: used to, one at a time, hold the names of the distribution
!              functions being read in.
      character(len=60) inputFile
!   numFiles: the number of files containing distribution function.
      integer numFiles
      real*8 average,energy
      real*8 p_bar(3)

      namelist /readDistValues/ mass,q,type_integration,numFiles,xstart,step

!  1.  read in the values from the namelist
      open(unit=9,file='readDistValues.dat',status='old')
      read(9,nml=readDistValues)
      close(9)

!  2.  open a file where the output data is to be printed
      open(unit=15,file='properties.out',status='unknown')
      ! hossna
      ! the output is vx,vy,vz which are devided (normalized with vth=sqrt(2.0 temp/mass)) and 
        ! df whch is also normalized to sqrt(2.0temp/mass)*sqrt(2.0temp/mass)*sqrt(2.0temp/mass), he 
        ! did that to avoid big or small numbers in output
        ! for the moments which are bigger than zero (1- (v df d3v ) or 2-( vv df d3v )or 3- ( vvv df d3v )and ...) depends the 
        ! number of v is used in the integration we should multiply the normalization factor to the results. 
        ! for the zero moment we dont need any normalization because the n= intgral (df d3v) so the df is normalized to cubic vth and v is normalized to vth
        ! and they cancle each other.
        ! but for other moments like 1- we should multiply the results by vth
        !                             2- we should multiply the results by vth*vth
        !                            3- we should multiply the results by vth*vth*vth

      average=0.
!  3.  do calculations for each of the files, one by one.
      dnintegr=0.
      pintegr=0.
      do i=2,numFiles+1
!   3.1  read in txyzf, and tindex from the file.
         call readFile(i)
!   3.2  calculate the density
         call calculations(1,density)
       !  density=density
       ! if(i<=36)then
       !   average=average+density
       ! endif
!   3.3  calculate the average momentum/velocity (which one
!        depends on which integration scheme was used to 
!        calculate the particle trajectories.  If symplectic
!        integration was used, the average momentum is 
!        calculated, while if runge-kutta integration was
!        used, the average velocity is calculated.)
         call calculations(2,momentum)
!   Hossna
        
         Ti=5.7*q
         momentum=momentum*sqrt(2.0*Ti/mass)
!rm      call calculations(2,p_bar)
!rm      call calculations(2,momentum)
!   3.4  calculate the energy
!rm      call calculations(3,energy)

!   3.5  calculate the flux vector
         !!!!!!call calculations(4,momentum)

!   3.6  calculate the stress tensor
         !!!!!!call calculations(6,stressTensor)
         ! hossna
         !!stressTensor=stressTensor*(2.0*Ti/mass)

        !!!!! dnintegr=dnintegr+abs(step)*1.e-6*density
        !!!!! pintegr=pintegr+abs(step)*(stressTensor(1)+stressTensor(4) &
        !!!!!   +stressTensor(6))*1.e9/3.

!   4.0  print the results in the output file
!rm      write(15,102)xstart+(i-1.)*step,density,energy,(energy/density),&
!rm           p_bar(1),p_bar(2),p_bar(3)
!  multiply stressTensor by 2.e-12 to include the electron contribution
!  and have the pressure in nPa.
!        write(15,103)xstart+(i-1.)*step,1.e-6*density,momentum, &
!          1.e+09*stressTensor,dnintegr,pintegr
        !!! write(15,103)xstart+(i-1.)*step,density,momentum, &
        !!!   stressTensor,dnintegr,pintegr
            write(15,103)xstart-(i-2.)*step,density,momentum
!   5.0  deallocate the txyzf, and tindex arrays
         deallocate(txyzf)
         deallocate(tindex)
      enddo

!rm print*,'average=',average
102 format(7es16.6)
103 format(99es14.6)
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
      !fileNum=fileNum+1
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
!rm         txyzf(4,i)=txyzf(4,i)/(sqrt(2.*mass*q)*sqrt(2.*mass*q)*&
!rm              sqrt(2.*mass*q))
            txyzf(4,i)=txyzf(4,i)/sqrt(2.*mass*q)**3
         enddo
      else
         do i=1,nVertices
            txyzf(1,i)=txyzf(1,i)*sqrt(2.*q/mass)
            txyzf(2,i)=txyzf(2,i)*sqrt(2.*q/mass)
            txyzf(3,i)=txyzf(3,i)*sqrt(2.*q/mass)
!rm         txyzf(4,i)=txyzf(4,i)/(sqrt(2.*q/mass)*sqrt(2.*q/mass)*&
!rm              sqrt(2.*q/mass))
            txyzf(4,i)=txyzf(4,i)/sqrt(2.*q/mass)**3
         enddo
      endif
endif
    
      return
      end

!=======================================================================
      subroutine calculations(type,calcValue)
      use VALUES
      implicit none

!  this subroutine calculates various quantities from the distribution
!  functions.


!  arguments
!   type specifies the quantity to be calculated.
!   calcValue: to contain the quantity calculated.
      integer type
      real*8 calcValue(*)

!  local variables
      integer i,j
      real*8 f(8),integral,density,pxAv,pyAv,pzAv

!  2.  calculate various quantities depending on the value of type
!   2.1  if type==1, calculate the density
      calcValue(1)=0.
      if(type==1)then
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(1)=calcValue(1)+integral
         enddo

!   2.2  if type==2, calculate the average velocity/momentum
      elseif(type==2)then  
         density=0.
         calcValue(1:3)=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            density=density+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(1,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(1)=calcValue(1)+integral
         enddo
         calcValue(1)=calcValue(1)/density
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(2,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(2)=calcValue(2)+integral
         enddo
         calcValue(2)=calcValue(2)/density
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(3,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(3)=calcValue(3)+integral
         enddo
         calcValue(3)=calcValue(3)/density
!   2.3  if type==3, calculate the energy
      elseif(type==3)then
         density=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            density=density+integral
         enddo
         calcValue(1)=0.
         if(type_integration==1)then
            do i=1,nCells
               do j=1,8
                  f(j)=0.5*(1./mass)*((txyzf(1,tindex(j,i))*txyzf(1,tindex(j,i)))+&
                       (txyzf(2,tindex(j,i))*txyzf(2,tindex(j,i)))+&
                       (txyzf(3,tindex(j,i))*txyzf(3,tindex(j,i))))*&
                       txyzf(4,tindex(j,i))
               enddo
               call integrateCell(i,f,integral)
               calcValue(1)=calcValue(1)+integral
            enddo
!            calcValue(1)=calcValue(1)/density
         else
            do i=1,nCells
               do j=1,8
                  f(j)=0.5*mass*((txyzf(1,tindex(j,i))*txyzf(1,tindex(j,i)))+&
                       (txyzf(2,tindex(j,i))*txyzf(2,tindex(j,i)))+&
                       (txyzf(3,tindex(j,i))*txyzf(3,tindex(j,i))))*&
                       txyzf(4,tindex(j,i))
               enddo
               call integrateCell(i,f,integral)
               calcValue(1)=calcValue(1)+integral
            enddo
!            calcValue(1)=calcValue(1)/density
         endif
!   2.4  if type==4, calculate the average flux
      elseif(type==4)then  
if(1 == 1) then
         density=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            density=density+integral
         enddo
else
  density=1.
endif

         calcValue(1:3)=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(1,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(1)=calcValue(1)+integral/density
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(2,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(2)=calcValue(2)+integral/density
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(3,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(3)=calcValue(3)+integral/density
         enddo
!   2.5  if type==5, calculate the total (including dynamic) stress
      elseif(type==5)then  
         calcValue(1:6)=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(1,tindex(j,i))*txyzf(1,tindex(j,i))/mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(1)=calcValue(1)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(1,tindex(j,i))*txyzf(2,tindex(j,i))/mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(2)=calcValue(2)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(1,tindex(j,i))*txyzf(3,tindex(j,i))/mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(3)=calcValue(3)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(2,tindex(j,i))*txyzf(2,tindex(j,i))/mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(4)=calcValue(4)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(2,tindex(j,i))*txyzf(3,tindex(j,i))/mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(5)=calcValue(5)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(3,tindex(j,i))*txyzf(3,tindex(j,i))/mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(6)=calcValue(6)+integral
         enddo
!   2.6  if type==6, calculate the stress tensor
      elseif(type==6)then  
         calcValue(1:6)=0.
!   2.6.1 first calculate the average momentum
         density=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            density=density+integral
         enddo

         pxAv=0.
         pyAv=0.
         pzAv=0.
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(1,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            pxAv=pxAv+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(2,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            pyAv=pyAv+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=txyzf(3,tindex(j,i))*txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            pzAv=pzAv+integral
         enddo
         pxAv=pxAv/density
         pyAv=pyAv/density
         pzAv=pyAv/density

         do i=1,nCells
            do j=1,8
               f(j)=(txyzf(1,tindex(j,i))-pxAv)*(txyzf(1,tindex(j,i))-pxAv)*mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(1)=calcValue(1)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=(txyzf(1,tindex(j,i))-pxAv)*(txyzf(2,tindex(j,i))-pyAv)*mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(2)=calcValue(2)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=(txyzf(1,tindex(j,i))-pxAv)*(txyzf(3,tindex(j,i))-pzAv)*mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(3)=calcValue(3)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=(txyzf(2,tindex(j,i))-pyAv)*(txyzf(2,tindex(j,i))-pyAv)*mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(4)=calcValue(4)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=(txyzf(2,tindex(j,i))-pyAv)*(txyzf(3,tindex(j,i))-pzAv)*mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(5)=calcValue(5)+integral
         enddo
         do i=1,nCells
            do j=1,8
               f(j)=(txyzf(3,tindex(j,i))-pzAv)*(txyzf(3,tindex(j,i))-pzAv)*mass &
                 *txyzf(4,tindex(j,i))
            enddo
            call integrateCell(i,f,integral)
            calcValue(6)=calcValue(6)+integral
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

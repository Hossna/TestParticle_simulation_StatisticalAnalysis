      program diagno
      implicit none
!  read in the output produced by readDist and produce various
!  diagnostics.

!  variables
      integer, parameter :: nout=11
      integer :: i,j,k
      real, allocatable :: tout(:,:),tav(:,:),fav(:)
      real surfCur
!  tout: field of values computed by readDist: coordinate, density, ...
!  tav: temporary array to hold a subset of tout over which to calculate
!       averages
!  fav: average values for various functions computed from tav
!  surfCur: Surface current density at the shock front

      real, parameter :: qe=1.6022e-19,mp=1.67262158e-27
!  qe: electron charge
!  mp: proton mass

      integer nStAv,nlines
      real step,delta
!  fName: name of file containing the output from readDist. Unually the
!         is something like 'properties.out'
!  nStAv: number of steps used to calculate the average.
!         for a periodic function, this is the number of (assumed
!         uniform) intervals
!  step:  step in one of the coordinates
!  delta: Length of the delta parameter in SI. In the code, all lengths
!         are normilized in terms of delta
!  nlines: number of lines in fName
      character(len=132) :: fName
      namelist/input/fName,nStAv,step,delta,nlines

!  procedures

!  computation

!  0.  Initialize variables

!  1.  Read in the output from readDist
      open(unit=3,file='diagno.dat',status='old')
      read(3,NML=input)

      allocate(tout(nlines,nout+1))
      allocate(fav(nout+3))

      open(unit=2,file=fName,status='old')
      do i=1,nlines
        read(2,*)tout(i,1:nout)
      enddo

!  2.  Compute integrals
      surfCur=0.
      do i=1,nlines
        surfCur=surfCur+delta*step*tout(i,3)*qe/mp
        tout(i,nout+1)=surfCur
      enddo

!  3.  Compute averages
      open(unit=7,file='diagno.out',status='unknown')
      if(nStAv > nLines) then
        write(6,*)'nStAv=',nStAv,'  nLines=',nLines
        write(6,*)'program will stop'
        stop
      endif
      do i=nStAv,nlines
        fav=0. !vector
        do j=i-nStAv+1,i
          do k=1,nout+1
            fav(k)=fav(k)+tout(j,k)
          enddo
        enddo
        fav=fav/nStAv !vector
!  total perpendicular stress
        fav(nout+2)=(fav(6)+fav(9))/2.
!  total stress, including paralellel stress
        fav(nout+3)=(fav(6)+fav(9)+fav(11))/3.
        write(7,102)(fav(k),k=1,nout+3)
 102    format(25es15.6)
      enddo
      close(unit=7)

!  4.  End the program
      stop
      end
!=======================================================================
      real function lagIntrpl(tx,ty,nt,n,x)
!  use Lagrange interpolation to calculate a function from tabulated
!  values. The tabulated values are in tx, ty. Each contains nt values.
!  The interpolation is perfoemed using n tabulated points on either
!  side (within bounds) of x.
!  n.b: the values in tx are assumed to be in ascending order

!  arguments
      integer n,nt
      real, intent(in) :: tx(nt),ty(nt),x

!  local variables
      integer i

!  procedures

!  computation

!  1.  find the interval containing x
      it=nt-1
      do i=1,nt-1
        if(x <= tx(i)) then
          it=i
          exit
        endif
      enddo

!  2.  find the extrema in the indices
      n1=max(1,it-n)
      n2=min(nt,it+n)

!  3.  calculate the inter(extra)polation with the Lagrange formula
     lagIntrpl=0.
     do i=n1,n2
       sum=1.
       do j=n1,n2
         if(i == j) cycle
         sum=sum*(x-tx(j))/(tx(i)-tx(j))
       enddo
       lagIntrpl=lagIntrpl+ty(i)*sum
     enddo
     return
     end

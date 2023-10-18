program readDist
! Copyritht 2009 by Richard Marchand,
! Department of Physics, University of Alberta.
! All rights reserved. This program or any part thereof may not be
! copied, modified or transferred to other parties without the explicit
! written consent of the author.
! The program may be used as is, and the author may not be held
! responsible for any loss or damage caused by its use.
!I used this program and made changes to simulate my own case study. 
!I put comment with my name "Hossna" everywhere inside the code.

**DESCRIPTION**:
this program reads in a set of distribution functions, and then calculates
various quantities from them.  The distribution functions must be called
distfunct*.out, where * increases from 000001 in steps of 1.  The x
position corresponding to the first distribution function
(distfunct000001.out), and the step by which the x position changes from
one distribution function to the next are specified in the file
readDistValues.dat.   
This program  prints the results in the file properties.out, where the 
first column is the x positions, the second column is the density, the 
third column is the energy, the fourth column is the energy divided by 
the density, and the last three columns are the average momentum or 
velocity.


**INPUT FILES**:
distfunct*.out: the distribution functions produced by the program
		densityOct, with denisty_type==2.  These distribution
		functions must start with distfunct000001.out, and increase
		in steps of 1.


**OUTPUT FILES**:
properties.out: this contains the quantities calculated by integrating
		the distribution function.
		The first column is the x positions, the second column is 
		the density, the third column is the energy, the fourth 
		column is the energy divided by the density, and the last 
		three columns are the average momentum or velocity.

**PROGRAM FILES**:
readDist.f90
values.mdl
readDistValues.dat 

Hossna: readDistNew.f90 (I worked with it)
          the output of tk.f90 code which is used here is: distfunct000001.tecplot and
          we should change the name to distfunct000001.out
        ! the output is vx,vy,vz which are devided (normalized with vth=sqrt(2.0 temp/mass)) and 
        ! df whch is also normalized to sqrt(2.0temp/mass)*sqrt(2.0temp/mass)*sqrt(2.0temp/mass), he 
        ! did that to avoid big or small numbers in output
        ! for the moments which are bigger than zero (1- (v df d3v ) or 2-( vv df d3v )or 3- ( vvv df d3v )and ...) depends the 
        ! number of v is used in the integration we should multiply the normalization factor to the results. 
        ! for the zero moment we dont need any normalization because the n= intgral (df d3v) so the df is normalized to cubic vth and v is normalized to vth
        ! and they cancle each other.
        but for other moments like 1- we should multiply the results by vth
                                   2- we should multiply the results by vth*vth
                                   3- we should multiply the results by vth*vth*vth

OBJS=	readDistNew.o 

#FLAGS= -r8 -O0 -g -C
#FLAGS=	-O5 -q64 -b64 
#FLAGS=	-r8 -O3
#F90=	ifort

FLAGS=	-O0 -fdefault-double-8 -fdefault-real-8 -fcheck=all -fbounds-check
FLAGS=	-O3 -fdefault-double-8 -fdefault-real-8
F90=	gfortran

readDistNew:	$(OBJS)
	$(F90) -o readDistNew $(FLAGS) $(OBJS)

readDistNew.o:	readDistNew.f90
	$(F90) -c $(FLAGS) readDistNew.f90

clean:
	rm -f $(OBJS) readDistNew

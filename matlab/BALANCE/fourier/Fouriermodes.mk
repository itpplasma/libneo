#FC       = f95f
#FC       = f95-lah
FC       = gfortran

#OPTS= -M OBJS --chk a,e,s,u,x --trace --trap -g
#OPTS= -M OBJS -O
OPTS= -J OBJS -O

OBJS =  OBJS/field_divB0.o \
	OBJS/bdivfree.o \
	OBJS/spline5_RZ.o \
	OBJS/rk4dc.o \
	OBJS/mag.o \
	OBJS/rhs_flt.o \
	OBJS/fouriermodes.o \

fouriermodes.x: $(OBJS) Fouriermodes.mk
	$(FC) $(OPTS) -o fouriermodes.x $(OBJS)
OBJS/field_divB0.o: field_divB0.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c field_divB0.f90
	mv field_divB0.o OBJS/
OBJS/bdivfree.o: bdivfree.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c bdivfree.f90
	mv bdivfree.o OBJS/
OBJS/spline5_RZ.o: spline5_RZ.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c spline5_RZ.f90
	mv spline5_RZ.o OBJS/
OBJS/rk4dc.o: rk4dc.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c rk4dc.f90
	mv rk4dc.o OBJS/
OBJS/mag.o: mag.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c mag.f90
	mv mag.o OBJS/
OBJS/rhs_flt.o: rhs_flt.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c rhs_flt.f90
	mv rhs_flt.o OBJS/
OBJS/fouriermodes.o: fouriermodes.f90 Fouriermodes.mk
	$(FC) $(OPTS) -c fouriermodes.f90
	mv fouriermodes.o OBJS/
#OBJS/.o: .f90 Fouriermodes.mk
#	$(FC) $(OPTS) -c .f90
#	mv .o OBJS/

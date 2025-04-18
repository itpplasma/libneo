#FC       = f95f
#FC       = f95-lah
FC = gfortran

#OPTS= -J OBJS --chk a,e,s,u,x --trace --trap -g
#OPTS= -M OBJS -O
#OPTS= -std=gnu -M OBJS -O3
#OPTS= -M OBJS -Waliasing -Wampersand  -Wline-truncation  -Wnonstd-intrinsics  -Wsurprising -Wno-tabs  -Wunderflow #-Wall# -Wunused-parameter -Wconversion -Wimplicit-interface -Wcharacter-truncation
#OPTS=-Wunused
#OPTS = -J OBJS -O -fPIC
OPTS = -J OBJS -O -fPIC -ffpe-trap=invalid,zero,overflow

OBJS =  OBJS/efit_to_boozer_mod.o \
	OBJS/odeint_allroutines.o \
	OBJS/spl_three_to_five_mod.o \
	OBJS/spline5_RZ.o \
	OBJS/field_divB0.o \
	OBJS/bdivfree_coul.o \
	OBJS/plag_coeff.o \
	OBJS/binsrc.o \
	OBJS/rhs_converter.o \
	OBJS/field_line_integration_for_converter.o \
	OBJS/hamada_to_boozer.o

hamada_to_boozer.x: $(OBJS) Hamada_to_Boozer.mk
	$(FC) $(OPTS) -o hamada_to_boozer.x $(OBJS) -llapack
OBJS/efit_to_boozer_mod.o: SRC/efit_to_boozer_mod.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/efit_to_boozer_mod.f90
	mv efit_to_boozer_mod.o OBJS/
OBJS/field_divB0.o: SRC/field_divB0.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/field_divB0.f90
	mv field_divB0.o OBJS/
OBJS/bdivfree_coul.o: SRC/bdivfree_coul.f90 Hamada_to_Boozer.mk SRC/field_divB0.f90
	$(FC) $(OPTS) -c SRC/bdivfree_coul.f90
	mv bdivfree_coul.o OBJS/
OBJS/spline5_RZ.o: SRC/spline5_RZ.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/spline5_RZ.f90
	mv spline5_RZ.o OBJS/
OBJS/spl_three_to_five_mod.o: SRC/spl_three_to_five_mod.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/spl_three_to_five_mod.f90
	mv spl_three_to_five_mod.o OBJS/
OBJS/odeint_allroutines.o: SRC/odeint_allroutines.f Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/odeint_allroutines.f
	mv odeint_allroutines.o OBJS/
OBJS/plag_coeff.o: SRC/plag_coeff.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/plag_coeff.f90
	mv plag_coeff.o OBJS/
OBJS/binsrc.o: SRC/binsrc.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/binsrc.f90
	mv binsrc.o OBJS/
OBJS/rhs_converter.o: SRC/rhs_converter.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/rhs_converter.f90
	mv rhs_converter.o OBJS/
OBJS/field_line_integration_for_converter.o: SRC/field_line_integration_for_converter.f90 Hamada_to_Boozer.mk
	$(FC) $(OPTS) -c SRC/field_line_integration_for_converter.f90
	mv field_line_integration_for_converter.o OBJS/
OBJS/hamada_to_boozer.o: SRC/hamada_to_boozer.f90 Hamada_to_Boozer.mk SRC/efit_to_boozer_mod.f90
	$(FC) $(OPTS) -c SRC/hamada_to_boozer.f90
	mv hamada_to_boozer.o OBJS/
#OBJS/.o: SRC/.f90 Hamada_to_Boozer.mk
#	$(FC) $(OPTS) -c SRC/.f90
#	mv .o OBJS/

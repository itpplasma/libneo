OBJS =  OBJS/vvn_asdex.o \
	OBJS/kisslinger_asdex.o
#COMP = f95-lah
COMP = gfortran
#OPTS = -O -M OBJS
OPTS = -O -J OBJS
%OPTS = --chk a,e,s,u,x -M OBJS
kisslinger_asdex.x: $(OBJS)
	$(COMP) $(OPTS) -o kisslinger_asdex.x $(OBJS)
OBJS/vvn_asdex.o: vvn_asdex.f
	$(COMP) $(OPTS) -c vvn_asdex.f
	mv vvn_asdex.o OBJS
OBJS/kisslinger_asdex.o: kisslinger_asdex.f90
	$(COMP) $(OPTS) -c kisslinger_asdex.f90
	mv kisslinger_asdex.o OBJS

FC = gfortran
LIB_HOME = ../../SIMPLE/build
LIB_NAME = simple
SOURCES = pyvmec.f90 ../../SIMPLE/SRC/magfie.f90 \
          ../../SIMPLE/SRC/spline_vmec_data.f90

lib: pyvmec.pyf pyvmec.f90
	f2py -c $^ \
	    --f90flags="-g -fbacktrace -I$(LIB_HOME) -fPIC" \
	    --fcompiler=$(FC) \
	    -L$(LIB_HOME) -l$(LIB_NAME)

pyvmec.pyf: $(SOURCES)
	f2py -m pyvmec -h $@ $^ --overwrite-signature

all: lib

clean:
	rm -rf *.pyf pyvmec*.so

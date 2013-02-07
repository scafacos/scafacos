#export PKG_CONFIG_PATH=${HOME}/workspace/scafacos/install/lib/pkg-config

all: scafacos_verify scafacos_verify_f

FCSLIB = $(shell pkg-config --libs scafacos-fcs)
FCSINC = $(shell pkg-config --cflags scafacos-fcs)


%.o: %.c
	mpicc -g -c $(FCSINC) -o $@ $^

%.fo: %.f90
	mpif90 -g -cpp -c $(FCSINC) -o $@ $^

scafacos_verify: setup.o scafacos_verify.o
	mpicxx -o $@ $^ $(FCSLIB) 

scafacos_verify_f: setup.o scafacos_verify.fo
	mpif90 -ffree-line-length-none -o $@ $^ $(FCSLIB) 

clean:
	rm -rf  *.o scafacos_verify
	rm -rf *.fo scafacos_verify_f

test: all
	./verify.sh

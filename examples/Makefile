FC=mpif90
FC_FLAGS=-O3 -Wall
FCSLIB = `pkg-config --libs scafacos`
FCSINC = `pkg-config --cflags scafacos`
fmm_test:
# F90 due to the need for proprocessing the files
	${FC} ${FCSINC} ${FC_FLAGS} -c fmm_test.F90
# linking the program to the library 
	${FC} -o fmm_test fmm_test.o ${FCSLIB}
clean:
	rm fmm_test
	rm fmm_test.o

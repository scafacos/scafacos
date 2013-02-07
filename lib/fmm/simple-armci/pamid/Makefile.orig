#CC      := mpicc
#COPT    := -g -O0 -std=c99 -Wall
#CC      := mpixlc_r
#COPT    := -g -O3 -qstrict -qlanglvl=extc99
CC      := /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc64-bgq-linux-gcc
COPT    := -g -O0 -std=c99 -Wall

SPI_DIR  = /bgsys/drivers/ppcfloor/spi
#PAMI_DIR = /bgsys/drivers/ppcfloor/comm/sys
PAMI_DIR = /home/xajhammond/PAMI/patch-9nov2011/sys
MPI_DIR  = /bgsys/drivers/ppcfloor/comm/gcc
INC      = -I$(MPI_DIR)/include -I$(PAMI_DIR)/include -I$(SPI_DIR)/include
LIB      = -L$(MPI_DIR)/lib     -L$(PAMI_DIR)/lib     -L$(SPI_DIR)/lib     -lmpich -lopa -lmpl -lrt -ldl -lpami -lSPI -lSPI_cnk -lpthread -lrt

LD      = $(CC)
CFLAGS  = $(COPT) $(INC) -DDEBUG_FUNCTION_ENTER_EXIT -I. -I../generic
LDFLAGS = $(COPT)

TEST_BINARIES := $(patsubst %.c,%.x,$(wildcard test*.c))
TEST_OBJECTS  := $(patsubst %.c,%.o,$(wildcard test*.c))
A1D_OBJECTS   := $(patsubst %.c,%.o,$(wildcard a1d*.c))
A1D_HEADERS   := $(patsubst %.c,%.h,$(wildcard a1d*.c))
OBJECTS       := $(A1D_OBJECTS)
HEADERS       := $(A1D_HEADERS)

all: $(TEST_BINARIES)

%.x: %.o liba1d.a
	$(LD) $(LDFLAGS) $< liba1d.a $(LIB) -o $@

liba1d.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS) $(TEST_OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TEST_BINARIES) liba1d.a


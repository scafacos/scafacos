DEVICE  = NONE

CC      = mpicc
COPT    = -g -O0 -std=c99 -Wall

ifeq ($(TARGET),BGP)
CC      = mpicc
COPT    = -g -O0 -std=c99 -Wall -D__bgp__
#CC      = mpixlc_r
#COPT    = -g -O3 -qstrict -qlanglvl=extc99 -D__bgp__
INC     = -I. -DDCMFD
DEVICE  = dcmfd
endif

ifeq ($(TARGET),BGQ)
CC      = mpicc
COPT    = -g -O0 -std=c99 -Wall -D__bgq__
#CC      = mpixlc_r
#COPT    = -g -O3 -qstrict -qlanglvl=extc99 -D__bgq__
INC     = -I. -DPAMID
DEVICE  = pamid
endif

ifeq ($(TARGET),XE)
CC      = cc
COPT    = -g -O0 -std=c99
INC     = -I. -DDMAPPD
DEVICE  = dmappd
endif

ifeq ($(DEVICE),NONE) 
@ You need to set TARGET={BGP,BGQ,XE}
endif

LIBS    = libparmci.a $(DEVICE)/liba1d.a generic/liba1g.a

LD      = $(CC)
CFLAGS  = $(COPT) $(INC) -I./generic
LDFLAGS = $(COPT)

TEST_BINARIES := $(patsubst %.c,%.x,$(wildcard test*.c))
TEST_OBJECTS  := $(patsubst %.c,%.o,$(wildcard test*.c))
OBJECTS       := parmci.o
HEADERS       := parmci.h

all: liba1g.a liba1d.a libparmci.a $(TEST_BINARIES)

%.x: %.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

libparmci.a: parmci.o
	$(AR) $(ARFLAGS) $@ $<

liba1d.a:
	make -C $(DEVICE)

liba1g.a:
	make CC=$(CC) -C generic

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS) $(TEST_OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TEST_BINARIES) libparmci.a
	make -C $(DEVICE) realclean
	make -C generic realclean


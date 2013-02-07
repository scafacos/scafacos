# Environment variable A1_INSTALL should point to the install directory of A1
# Environment variable ARMCI_INSTALL should point to the install directory of ARMCI
# Environment variable GA_INSTALL should point to the install directory of GA

CC      = mpixlc_r 
COPT    = -g -O0

INCLUDE = -I$(GA_INSTALL)/include -I$(A1_INSTALL)/include 
LIBS    = -I$(GA_INSTALL)/lib -lglobal -ltcgmsg-mpi -lma -L$(A1_INSTALL)/lib -la1 

LD      = $(CC)
CFLAGS  = $(COPT) $(INCLUDE)
LDFLAGS = $(COPT)

TESTS := $(patsubst %.c,%.x,$(wildcard *.c))
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))

all: $(TESTS)

%.x: %.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TESTS)


# Environment variable A1_INSTALL should point to the install directory of A1 

CC      = mpixlc_r 
COPT    = -g -O0

INCLUDE = -I. -I$(A1_INSTALL)/include 
LIBS    = -L$(A1_INSTALL)/lib -la1 

LD      = $(CC)
CFLAGS  = $(COPT) $(INCLUDE)
LDFLAGS = $(COPT)

TESTS := $(patsubst %.c,%.x,$(wildcard *.c))
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))

all: $(TESTS)
	mv perf.x perf/.
	mv msgcheck.x msgcheck/.
	mv test-ibm.x test-ibm/.
	mv test2.x test2/.
	mv test.x test/.
	mv perf_nb.x perf_nb/.
	mv perf_aggr.x perf_aggr/.

%.x: %.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TESTS)


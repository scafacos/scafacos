FC      = ftn 
CC      = cc 
COPT    = -g -O2
FOPT    = -g -O2

# CSCS paths
DMAPP_PATH=/opt/cray/dmapp/2.2-1.0301.2349.6.6.gem
GNIH_PATH=/opt/cray/gni-headers/2.0-1.0301.2302.4.1.gem
PMI_PATH=/opt/cray/pmi/1.0-1.0000.7901.22.1.gem

INCLUDE = -I$(DMAPP_PATH)/include -I$(GNIH_PATH)/include -I$(PMI_PATH)/include
LIBS    = -L$(DMAPP_PATH)/lib64   -L$(GNIH_PATH)/lib64   -L$(PMI_PATH)/lib64 -ldmapp -lpmi

CFLAGS  = $(COPT) $(INCLUDE)

LD      = $(CC)
LDFLAGS = $(COPT)

TESTS := $(patsubst %.c,%.x,$(wildcard *.c))
OBJ   := $(patsubst %.c,%.o,$(wildcard *.c))
ASM   := $(patsubst %.c,%.s,$(wildcard *.c))

all: tests

tests: $(TESTS)

%.x: %.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.s: %.c
	$(CC) $(CFLAGS) -S $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ)

realclean: clean
	$(RM) $(RMFLAGS) $(TESTS)
	$(RM) $(RMFLAGS) $(ASM)


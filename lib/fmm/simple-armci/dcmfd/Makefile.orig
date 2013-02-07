CC      = mpicc
COPT    = -g -O0 -std=gnu99
#CC      = mpixlc_r
#COPT    = -g -O3 -qstrict -qlanglvl=extc99

INC     = -I. -I../generic -I/bgsys/drivers/ppcfloor/arch/include/

LD      = $(CC)
CFLAGS  = $(COPT) $(INC) -DDEBUG_FUNCTION_ENTER_EXIT
LDFLAGS = $(COPT)

TEST_BINARIES := $(patsubst %.c,%.x,$(wildcard test*.c))
TEST_OBJECTS  := $(patsubst %.c,%.o,$(wildcard test*.c))
A1D_OBJECTS   := $(patsubst %.c,%.o,$(wildcard a1d*.c))
A1D_HEADERS   := $(patsubst %.c,%.h,$(wildcard a1d*.c))
OBJECTS       := $(A1D_OBJECTS)
HEADERS       := $(A1D_HEADERS)

all: $(TEST_BINARIES)

%.x: %.o liba1d.a
	$(LD) $(LDFLAGS) $< liba1d.a -o $@

liba1d.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS) $(TEST_OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TEST_BINARIES) liba1d.a


CC      = mpixlc_r 
COPT    = -g -O0

INCLUDE = -I. -I/bgsys/drivers/ppcfloor/arch/include/
LIBS    = libbench.a

LD      = $(CC)
CFLAGS  = $(COPT) $(INCLUDE)
LDFLAGS = $(COPT)

GET_TESTS := $(patsubst %.c,%.x,$(wildcard get*.c))
PUT_TESTS := $(patsubst %.c,%.x,$(wildcard put*.c))
SEND_TESTS := $(patsubst %.c,%.x,$(wildcard send*.c))
BGP_TESTS := $(patsubst %.c,%.x,$(wildcard bgp*.c))
COLL_TESTS := $(patsubst %.c,%.x,$(wildcard many*.c))

TESTS := $(GET_TESTS) $(PUT_TESTS) $(SEND_TESTS) $(BGP_TESTS) $(COLL_TESTS)
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))

all: libbench.a $(TESTS)

libbench.a: bench_impl.o
	$(CC) $(CFLAGS) -c $< -o $@
	$(AR) $(ARFLAGS) libbench.a bench_impl.o

%.x: %.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TESTS) libbench.a


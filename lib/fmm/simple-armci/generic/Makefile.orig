CC      := mpicc
COPT    := -g -O0 -std=c99

INC     = -I.

LD      = $(CC)
CFLAGS  = $(COPT) $(INC)
LDFLAGS = $(COPT)

TEST_BINARIES := $(patsubst %.c,%.x,$(wildcard test*.c))
TEST_OBJECTS  := $(patsubst %.c,%.o,$(wildcard test*.c))
A1G_OBJECTS   := $(patsubst %.c,%.o,$(wildcard a1d*.c))
A1G_HEADERS   := $(patsubst %.c,%.h,$(wildcard a1d*.c))
OBJECTS       := $(A1G_OBJECTS)
HEADERS       := $(A1G_HEADERS)

all: liba1g.a $(TEST_BINARIES)

%.x: %.o liba1g.a
	$(LD) $(LDFLAGS) $< liba1g.a -o $@

liba1g.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJECTS) $(TEST_OBJECTS)

realclean: clean
	$(RM) $(RMFLAGS) $(TEST_BINARIES) liba1g.a


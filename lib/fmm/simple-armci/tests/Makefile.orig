CC      = mpicc
COPT    = -g -O2 -Wall -std=c99
#CC      = mpixlc_r
#COPT    = -g -O3 -qstrict -qlanglvl=extc99

INC     = -I../
LIBS    = -L../ -larmci-bgp

LD      = $(CC)
CFLAGS  = $(COPT) $(INC)
LDFLAGS = $(COPT)

all: ARMCI_Hello.x ARMCI_Get.x

#ARMCI_Get_bw.x ARMCI_Get_latency.x ARMCI_Put_bw.x ARMCI_Put_latency.x

%.x: %.o
	$(LD) $(LDFLAGS) $< $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.x


.PHONE = all clean

CC = gcc
CXX = gcc
CFLAGS = -O
CPPFLAGS = -O
LDFLAGS = -lm -lstdc++
OBJECT_FILES = gendata_test.o gendata_util.o gendata_hammersley.o hammersley.o \
               gendata_grid.o gendata_plummer.o

all: gendata

clean:
	rm -f gendata
	rm -f $(OBJECT_FILES)

gendata: $(OBJECT_FILES)
	$(CC) -o $@ $(OBJECT_FILES) $(LDFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) -o $@ $<


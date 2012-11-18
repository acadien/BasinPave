CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(notdir $(CPP_FILES:.cpp=.o))
INC     = -I/usr/include/c++/4.6 
LIBS 	= -lgsl -lgslcblas -lm -lstdc++ -llapack
CC	= gcc -O3

basinPave: $(OBJ_FILES) nrutil.o
	$(CC) -g $(INC) -o $@ $^ $(LIBS)

nrutil.o: nrutil.c
	$(CC) -g -c nrutil.c

%.o: %.cpp
	$(CC) -g $(LIBS) -c -o $@ $<

clean:
	rm *.o basinPave

all: 	basinPave
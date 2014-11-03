UNAME_M = $(shell uname -m)
UNAME_S = $(shell uname -s)

# -g is for debugging symbols
# CFALL = $(CFLAGS) $(TARGET_ARCH) -g 
CFALL = $(CFLAGS) $(TARGET_ARCH)  
ifeq ($(UNAME_M), ppc64)
	CFARCHDEP = -q64 -qtune=pwr7 -qarch=pwr7 -qhot
	MPCC = mpcc
	INC = 
	LIB = 
endif
ifeq ($(UNAME_S), Darwin)
	CFARCHDEP = -Wall -std=c99 -g
	MPCC = mpicc
	INC = -I/usr/local/include
	LIB = -L/usr/local/lib
else
	CFARCHDEP = -Wall -std=c99 -g
	MPCC = mpicc
	INC = -I/usr/include/suitesparse 
	LIB =
endif

CF = $(CFALL) $(CFARCHDEP)
CS = -lcxsparse
OBJ = matops.o brain.o nvu.o adjacency.o 
EXE = testmatops testmat2 simulate testbrain

all: tags $(OBJ) $(EXE)

tags: *.h *.c
	ctags *.h *.c

brain.o: brain.c brain.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c brain.c 

adjacency.o: adjacency.c brain.h Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c adjacency.c 

nvu.o: nvu.h nvu.c Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c nvu.c

matops.o: matops.h matops.c Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -c matops.c

testmatops: testmatops.c matops.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o testmatops testmatops.c matops.o $(CS) -lm 

testmat2: testmat2.c matops.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o testmat2 testmat2.c matops.o $(CS) -lm

testbrain: testbrain.c matops.o brain.o adjacency.o nvu.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o testbrain testbrain.c $(OBJ) $(CS) -lm

simulate: simulate.c matops.o brain.o nvu.o Makefile
	$(MPCC) $(CF) $(INC) $(LIB) -o simulate simulate.c $(OBJ) $(CS) -lm

clean:
	rm -r tags $(OBJ) $(EXE) *.dSYM

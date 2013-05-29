MACHINE = $(shell uname -m)

CFALL = $(CFLAGS) $(TARGET_ARCH) -g 
ifeq ($(MACHINE),ppc64)
	CFARCHDEP = -q64 -qtune=pwr7 -qarch=pwr7 -qhot
	MPCC = mpcc
	INC = 
else
	CFARCHDEP = -Wall -std=c99
	MPCC = mpicc
	INC = -I/usr/include/suitesparse
endif

CF = $(CFALL) $(CFARCHDEP)
CS = -lcxsparse
OBJ = matops.o brain.o nvu.o adjacency.o
EXE = testmatops testmat2 simulate testbrain

all: tags $(OBJ) $(EXE)

tags: *.h *.c
	ctags *.h *.c

brain.o: brain.c brain.h Makefile
	$(MPCC) $(CF) $(INC) -c brain.c 

adjacency.o: adjacency.c brain.h Makefile
	$(MPCC) $(CF) $(INC) -c adjacency.c 

#constants.o: constants.c brain.h Makefile
#	$(MPCC) $(CF) $(INC) -c constants.c

nvu.o: nvu.h nvu.c Makefile
	$(MPCC) $(CF) $(INC) -c nvu.c

matops.o: matops.h matops.c Makefile
	$(MPCC) $(CF) $(INC) -c matops.c

testmatops: testmatops.c matops.o Makefile
	$(MPCC) $(CF) $(INC) -o testmatops testmatops.c matops.o $(CS) -lm 

testmat2: testmat2.c matops.o Makefile
	$(MPCC) $(CF) $(INC) -o testmat2 testmat2.c matops.o $(CS) -lm

testbrain: testbrain.c matops.o brain.o adjacency.o Makefile
	$(MPCC) $(CF) $(INC) -o testbrain testbrain.c $(OBJ) $(CS) -lm

simulate: simulate.c matops.o brain.o nvu.o Makefile
	$(MPCC) $(CF) $(INC) -o simulate simulate.c $(OBJ) $(CS) -lm

clean:
	rm tags $(OBJ) $(EXE)

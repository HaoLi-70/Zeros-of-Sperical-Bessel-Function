# Makefile MCMC inversion
#-------------------------------------------
CC = gcc
ARRCS = ar rcs
CFLAGS = -c

all: Run libmine.a clean

OBJ1=ALLOCATION.o SPECIAL_FUNCTIONS.o TIME_PRINT.o\
       
OBJ2=Compute_Zeros.o	

ALLOCATION.o:ALLOCATION.c
	$(CC) $(CFLAGS) ALLOCATION.c

SPECIAL_FUNCTIONS.o:SPECIAL_FUNCTIONS.c 
	$(CC) $(CFLAGS) SPECIAL_FUNCTIONS.c

TIME_PRINT.o:TIME_PRINT.c
	$(CC) $(CFLAGS) TIME_PRINT.c

Compute_Zeros.o:Compute_Zeros.c
	$(CC) $(CFLAGS) Compute_Zeros.c

libmine.a:$(OBJ1)
	$(ARRCS) libmine.a $(OBJ1)

Run: $(OBJ2) libmine.a
	$(CC) $(OBJ2) -o Run -L. -lm -lmine

.PHONY.:clean

clean:
	rm -f *.o *.a



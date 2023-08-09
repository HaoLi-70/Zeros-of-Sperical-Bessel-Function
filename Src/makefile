# Makefile MCMC inversion
#-------------------------------------------


CMP = gcc
CMPS = gnu

ifeq ($(CMPS), gnu)
ifeq ($(debug), y)
	CFLAGS := -c -Wall -Wshadow -Wextra -Warray-bounds -pedantic
else
ifeq ($(hard), y)
	CFLAGS := -c -O3
else
	CFLAGS := -c
endif
endif
endif

all: Zero

SRC := $(wildcard *.c)
OBJS := $(SRC:.c=.o)
#OBJ := $(patsubst %.c,%.o,$(wildcard *.c))

%.o: %.c
	$(CMP) $(CFLAGS) -o $@ $< 

Zero: $(OBJS)
	$(CMP) -o Zero $(OBJS) -L. -lm
	mv Zero ../

.PHONY:clean

clean:
	find . -maxdepth 2 -name "*.o" -delete 
	find . -maxdepth 2 -name "Zero" -delete 
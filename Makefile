MCBSPLIB=../libs/mcbsp/lib/libmcbsp1.1.0.a
MCBSPINCLUDE=../libs/mcbsp/

CC=gcc 
CFLAGS=-O3 -std=c99 -I $(MCBSPINCLUDE)
LFLAGS=$(MCBSPLIB) -pthread -lm -lrt -lhwloc
OBJ= mbspinprod.o mbsp-discover.o mbsputil.o

all: mbspinprod 


mbspinprod: $(OBJ)
	$(CC) $(CFLAGS) -o mbspinprod $(OBJ) $(LFLAGS)

mbspinprod.o: mbspinprod.c 
	$(CC) $(CFLAGS) -c mbspinprod.c

mbsp-discover.o: mbsp-discover.c mbsp-discover.h
	$(CC) $(CFLAGS) -c mbsp-discover.c

mbsputil.o: mbsputil.c  mbsputil.h
	$(CC) $(CFLAGS) -c mbsputil.c

clean:
	rm -f *.o mbspinprod 

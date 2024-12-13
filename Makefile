all: 	edf_hull
CFLAGS = -g -O0 -std=c89 -Wpedantic -D_GNU_SOURCE

ts_lib.o:	ts_lib.c ts_lib.h Makefile
	$(CC) -c $(CFLAGS)  ts_lib.c -o ts_lib.o

edf_hull.o: edf_hull.c edf_hull.h Makefile
	$(CC) -c $(CFLAGS)  edf_hull.c -o edf_hull.o

edf_hull: edf_hull.o ts_lib.o main.c main.h  Makefile
	$(CC) -g -O0  main.c edf_hull.o ts_lib.o -lglpk -lm -o edf_hull

clean:
	rm -rf *.o *~ edf_hull edf_debug_qhull.txt

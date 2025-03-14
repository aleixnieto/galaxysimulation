CC = gcc
CFLAGS = -Wall -O3 -march=native -fopt-info-vec -ffast-math -funroll-loops
LDFLAGS = -lm $(shell pkg-config --libs x11)

all: galsim compare_gal_files

galsim: galsim.c graphics.c graphics.h
	$(CC) $(CFLAGS) galsim.c graphics.c -o galsim $(LDFLAGS)

compare_gal_files: compare_gal_files.c
	$(CC) $(CFLAGS) compare_gal_files.c -o compare_gal_files $(LDFLAGS)

clean:
	rm -f galsim compare_gal_files

clean_results:
	rm -f result.gal

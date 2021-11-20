CC=gcc
CFLAGS = -O2 -lm -fopenmp -fPIC -Wall -m64

default: all

triangleCount_openMP: triangleCount_openMP.c
	$(CC) $(CFLAGS) -o $@ $^ mmio.c

.PHONY: clean

all: 	triangleCount_openMP

test: triangleCount_openMP
	echo "binary,graph,node_num,threads_num,times,triang_num" > openMP_runtimes.csv
	@printf "\n** Testing triangleCount_openMP\n\n"
	./triangleCount_openMP belgium_osm.mtx
	./triangleCount_openMP com-Youtube.mtx
	./triangleCount_openMP mycielskian13.mtx
	./triangleCount_openMP dblp-2010.mtx
	./triangleCount_openMP NACA0015.mtx

clean:
	rm -f triangleCount_openMP

CC=gcc
CFLAGS = -O2 -Wall -lm -m64 -pthread -fPIC

default: all

triangleCount_Pthreads: triangleCount_Pthreads.c
	$(CC) $(CFLAGS) -o $@ $^ mmio.c

.PHONY: clean

all: 	triangleCount_Pthreads

test: triangleCount_Pthreads
	echo "binary,graph,node_num,threads_num,times,triang_num" > pthreads_runtimes.csv
	@printf "\n** Testing triangleCount_Pthreads\n\n"
	./triangleCount_Pthreads belgium_osm.mtx
	./triangleCount_Pthreads com-Youtube.mtx
	./triangleCount_Pthreads mycielskian13.mtx
	./triangleCount_Pthreads dblp-2010.mtx
	./triangleCount_Pthreads NACA0015.mtx

clean:
	rm -f triangleCount_Pthreads

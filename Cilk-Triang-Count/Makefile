CC=gcc
CFLAGS = -O2 -fPIC -m64 -Wall  -lm -lcilkrts -fcilkplus

default: all

triangleCount_Cilk: triangleCount_Cilk.c
	$(CC) $(CFLAGS) -o $@ $^ mmio.c

.PHONY: clean

all: 	triangleCount_Cilk

test: triangleCount_Cilk
	echo "binary,graph,node_num,threads_num,times,triang_num" > cilk_runtimes.csv
	@printf "\n** Testing triangleCount_Cilk\n\n"
	./triangleCount_Cilk belgium_osm.mtx
	./triangleCount_Cilk com-Youtube.mtx
	./triangleCount_Cilk mycielskian13.mtx
	./triangleCount_Cilk dblp-2010.mtx
	./triangleCount_Cilk NACA0015.mtx

clean:
	rm -f triangleCount_Cilk

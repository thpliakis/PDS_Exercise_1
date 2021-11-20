CC=gcc
CFLAGS = -O3 -fPIC

default: all

triangleCount: triangleCount.c
	$(CC) $(CFLAGS) -o $@ $^ mmio.c

.PHONY: clean

all: triangleCount

test: triangleCount
	echo "binary,graph,node_num,threads_num,times,triang_num" > seq_runtimes.csv
	@printf "\n** Testing triangleCount\n\n"
	./triangleCount belgium_osm.mtx
	./triangleCount com-Youtube.mtx
	./triangleCount mycielskian13.mtx
	./triangleCount dblp-2010.mtx
	./triangleCount NACA0015.mtx
clean:
	rm -f triangleCount

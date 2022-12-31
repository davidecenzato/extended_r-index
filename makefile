# compilation flags
CCX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG -g -I /usr/local/include/ -L /usr/local/lib/
#CCX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -g -I /usr/local/include/ -L /usr/local/lib/
CCX=g++

# main executables 
EXECS = er-index er-index64 genpattern

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

er-index: main.cpp r_index.hpp rle_ebwt.hpp pred_ebwt.hpp sd_vector.hpp external/malloc_count/malloc_count.o 
	$(CCX) $(CCX_FLAGS) -o $@ main.cpp r_index.hpp rle_ebwt.hpp pred_ebwt.hpp sd_vector.hpp external/malloc_count/malloc_count.o -ldl -lsdsl -ldivsufsort -ldivsufsort64

er-index64: main.cpp r_index.hpp rle_ebwt.hpp pred_ebwt.hpp sd_vector.hpp external/malloc_count/malloc_count.o 
	$(CCX) $(CCX_FLAGS) -o $@ main.cpp r_index.hpp rle_ebwt.hpp pred_ebwt.hpp sd_vector.hpp external/malloc_count/malloc_count.o -ldl -lsdsl -ldivsufsort -ldivsufsort64 -DM64

genpattern: genpattern.cpp external/malloc_count/malloc_count.o 
	$(CCX) $(CCX_FLAGS) -o $@ genpattern.cpp external/malloc_count/malloc_count.o -ldl -lsdsl -ldivsufsort -ldivsufsort64

clean:
	rm -f $(EXECS)
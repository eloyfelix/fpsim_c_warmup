.PHONY: all clean
.SUFFIXES: .cpp .o .c .h

CFLAGS= -fPIC -std=gnu99 -Wall -Wextra -Wshadow

CFLAGS += -O3 -funroll-loops

CFLAGS_GCC =
CFLAGS_ICC = 

CFLAGS += -march=native #-DHAVE_AVX2_INSTRUCTIONS
# CFLAGS_GCC += -mavx2
# CFLAGS_ICC += -march=core-avx2

all: test_chembl sim_search.so 

HEADERS=./include/sim_search.h \
		./include/score_functions.h 
		
J_O=score_functions.o

SS_O=sim_search.o

%.o: ./src/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -Iinclude

sim_search.so: ./src/sim_search.c 
	$(CC) $(CFLAGS) -shared -o $@ $? -lc -Iinclude $(J_O)

test_chembl: ./src/test_chembl.c $(HEADERS) $(J_O) $(SS_O)
	$(CC) $(CFLAGS) -o $@ ./src/test_chembl.c -Iinclude  $(J_O) $(SS_O)

clean:
	rm -f test_chembl *.o *.so out.bin

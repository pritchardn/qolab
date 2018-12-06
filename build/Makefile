# the compiler: gcc for C program, define as g++ for C++
CC = icc

# compiler flags:
CFLAGS = -std=c99 -DMKL_ILP64 -O3 -I${MKLROOT}/include -fopenmp -Wall -Werror
LINKERS =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lnlopt -lm -ldl
# the build target executable:
TARGET = ../bin/qaoa.exe
LOC = ../src
SRCS = $(LOC)/main.c $(LOC)/qaoa.c $(LOC)/ub.c $(LOC)/globals.c
HEADERS = $(LOC)/main.h $(LOC)/qaoa.h $(LOC)/ub.h $(LOC)/globals.h

build: $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(HEADERS) $(LINKERS)

clean:
	rm -f *.o
	rm $(TARGET)

rebuild: clean build

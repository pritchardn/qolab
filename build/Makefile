# the compiler: gcc for C program, define as g++ for C++
CC = icc

# compiler flags:
CFLAGS = -std=c99 -DMKL_ILP64 -O3 -I${MKLROOT}/include -fopenmp -Wall -Werror
LINKERS =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lnlopt -lm -ldl
# the build target executable:
TARGET = ../bin/qaoa.exe
LOC = ../src
SRCS = $(LOC)/main.c $(LOC)/qaoa.c $(LOC)/ub.c $(LOC)/globals.c $(LOC)/uc.c $(LOC)/problem_code.c $(LOC)/state_evolve.c $(LOC)/reporting.c $(LOC)/matrix_expm.c $(LOC)/graph_utils.c
HEADERS = $(LOC)/qaoa.h $(LOC)/ub.h $(LOC)/globals.h $(LOC)/uc.h $(LOC)/problem_code.h $(LOC)/state_evolve.h $(LOC)/reporting.h $(LOC)/matrix_expm.h $(LOC)/graph_utils.h

build: $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(HEADERS) $(LINKERS)

clean:
	rm -f *.o
	rm $(TARGET)

rebuild: clean build

# the compiler: gcc for C program, define as g++ for C++
CC = icc

# compiler flags:
CFLAGS = -std=c99 -DMKL_LP64 -O3 -I${MKLROOT}/include -fopenmp -Wall -Werror
LINKERS =  -Wl, -mkl=parallel -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lnlopt
# the build target executable:
TARGET = ../bin/qaoa.exe
LOC = ../src
SRCS = $(LOC)/main.c $(LOC)/qaoa.c $(LOC)/ub.c $(LOC)/globals.c $(LOC)/uc.c $(LOC)/problem_code.c $(LOC)/state_evolve.c $(LOC)/reporting.c $(LOC)/matrix_expm.c $(LOC)/graph_utils.c $(LOC)/measurement.c $(LOC)/eigen_solve.c
HEADERS = $(LOC)/qaoa.h $(LOC)/ub.h $(LOC)/globals.h $(LOC)/uc.h $(LOC)/problem_code.h $(LOC)/state_evolve.h $(LOC)/reporting.h $(LOC)/matrix_expm.h $(LOC)/graph_utils.h $(LOC)/measurement.h $(LOC)/eigen_solve.h

build: $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(HEADERS) $(LINKERS)

clean:
	rm -f *.o
	rm $(TARGET)

rebuild: clean build

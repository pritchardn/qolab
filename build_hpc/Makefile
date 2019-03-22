# the compiler:
CC = cc

MAALI_NLOPT_HOME = /group/pawsey0309/npritchard/software/cle60up05/apps/PrgEnv-intel/6.0.4/intel/17.0.4.196/haswell/nlopt/2.5.0

# compiler flags: assumes nlopt is installed locally under $HOME/install
CFLAGS = -std=c99 -DMKL_LP64 -mkl=parallel -O3 -I$(MAALI_NLOPT_HOME)/include -Wall
LINKERS =  -liomp5 -lpthread -L$(MAALI_NLOPT_HOME)/lib64 -lnlopt

# the build target executable:
TARGET = ../bin/qaoa.exe
LOC = ../src
SRCS = $(LOC)/main.c $(LOC)/qaoa.c $(LOC)/ub.c $(LOC)/globals.c $(LOC)/uc.c $(LOC)/problem_code.c $(LOC)/state_evolve.c $(LOC)/reporting.c $(LOC)/matrix_expm.c $(LOC)/graph_utils.c $(LOC)/measurement.c $(LOC)/eigen_solve.c
HEADERS = $(LOC)/qaoa.h $(LOC)/ub.h $(LOC)/globals.h $(LOC)/uc.h $(LOC)/problem_code.h $(LOC)/state_evolve.h $(LOC)/reporting.h $(LOC)/matrix_expm.h $(LOC)/graph_utils.h $(LOC)/measurement.h $(LOC)/eigen_solve.h

build: $(SRCS)
	$(CC) $(CFLAGS) $(LINKERS) -o $(TARGET) $(SRCS) $(HEADERS) 

clean:
	rm -f *.o
	rm $(TARGET)

rebuild: clean build


ifdef CPLEX
        CPLEX_DIR = /opt/ibm/ILOG/CPLEX_Studio1271/cplex
        CPLEXFLAGS1 = -L$(CPLEX_DIR)/lib/x86-64_linux/static_pic -lcplex
        CPLEXFLAGS2 = -DCPLEX -I$(CPLEX_DIR)/include/ilcplex
        GLPKFLAGS =
else
        CPLEXFLAGS1 = 
        CPLEXFLAGS2 = 
        GLPKFLAGS = -lglpk
endif

HOMEBIN=.


CC=g++

#
# This can be left unchanged
#


LDFLAGS= -O3 -g $(CPLEXFLAGS1) $(GLPKFLAGS) -lpthread -lm -std=c++11 -Wall -Wextra -Wuninitialized -Wshadow
CFLAGS=  -O3 -g $(CPLEXFLAGS2) -std=c++11 -Wall -Wextra -Wuninitialized -Wshadow


OBJ =  	    multiobjective_solver.o problem_class.o simplex_class.o sydneys_class.o point_class.o SimplexStore.o
HEADER = 	multiobjective_solver.h problem_class.h simplex_class.h sydneys_class.h point_class.h SimplexStore.h

all: ${HOMEBIN}/MOS

clean:
	@rm *.o
	@[ -f ${HOMEBIN}/MOS ]

${HOMEBIN}/MOS: $(OBJ) Makefile $(HEADER)
	@echo Linking $(@F)
	@$(CC) -o ${HOMEBIN}/MOS $(OBJ) $(LDFLAGS)

%.o: %.c Makefile $(HEADER)
	@echo [${CC}] $<
	@$(CC) $(CFLAGS) -c $< 

#
# Set these two directories
#

HOMEBIN=${HOME}/Dropbox\ \(Edinboro\ University\)/Edinboro/Teaching/Student_Projects/Lesseski/solution_compare
CPLEX_DIR=$(HOME)/ILOG/CPLEX_Studio126/cplex

CC=g++

#
# This can be left unchanged
#

SOLSTRUCT=tree

LDFLAGS= -O3 -g -L$(CPLEX_DIR)/lib/x86-64_linux/static_pic -lcplex -lpthread -lm -std=c++11 -Wall -Wextra -Wuninitialized -Wshadow
CFLAGS=  -O3 -g -I$(CPLEX_DIR)/include/ilcplex -std=c++11 -Wall -Wextra -Wuninitialized -Wshadow

OBJ =  	    compare.o
HEADER = 	compare.h

all: ${HOMEBIN}/comp

clean:
	@rm *.o
	@[ -f ${HOMEBIN}/comp ]

${HOMEBIN}/comp: $(OBJ) Makefile $(HEADER)
	@echo Linking $(@F)
	@$(CC) -DSOL_$(SOLSTRUCT) -o ${HOMEBIN}/comp $(OBJ) $(LDFLAGS)

%.o: %.c Makefile $(HEADER)
	@echo [${CC}] $<
	@$(CC) -DSOL_$(SOLSTRUCT) $(CFLAGS) -c $< 

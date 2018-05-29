#
# Set these two directories
#

HOMEBIN=${HOME}/Dropbox\ \(Edinboro\ University\)/Edinboro/Teaching/Student_Projects/Lesseski/converter
CPLEX_DIR=$(HOME)/ILOG/CPLEX_Studio126/cplex

CC=g++

#
# This can be left unchanged
#

SOLSTRUCT=tree

LDFLAGS= -O3 -g -L$(CPLEX_DIR)/lib/x86-64_linux/static_pic -lcplex -lpthread -lm -std=c++11 -Wall -Wextra -Wuninitialized -Wshadow
CFLAGS=  -O3 -g -I$(CPLEX_DIR)/include/ilcplex -std=c++11 -Wall -Wextra -Wuninitialized -Wshadow

OBJ =  	    converter.o
HEADER = 	converter.h

all: ${HOMEBIN}/convert

clean:
	@rm *.o
	@[ -f ${HOMEBIN}/convert ]

${HOMEBIN}/convert: $(OBJ) Makefile $(HEADER)
	@echo Linking $(@F)
	@$(CC) -DSOL_$(SOLSTRUCT) -o ${HOMEBIN}/convert $(OBJ) $(LDFLAGS)

%.o: %.c Makefile $(HEADER)
	@echo [${CC}] $<
	@$(CC) -DSOL_$(SOLSTRUCT) $(CFLAGS) -c $< 

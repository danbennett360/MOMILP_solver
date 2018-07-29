#
# Set these two directories
#

HOMEBIN=${HOME}/Dropbox\ \(Edinboro\ University\)/Edinboro/Teaching/Student_Projects/Lesseski/Test_Instances/BelottiRan/generator/
#CPLEX_DIR=$(HOME)/ILOG/CPLEX_Studio126/cplex

CC=gcc

#
# This can be left unchanged
#

#SOLSTRUCT=tree

LDFLAGS= -O3 -g -lpthread -lm
CFLAGS=  -O3 -g 

OBJ = test_problems.o
HEADER = 

all: ${HOMEBIN}/generate_instances

clean:
	@rm *.o
	@[ -f ${HOMEBIN}/generate_instances ]

${HOMEBIN}/generate_instances: $(OBJ) Makefile $(HEADER)
	@echo Linking $(@F)
	@$(CC) -DSOL_$(SOLSTRUCT) -o ${HOMEBIN}/generate_instances $(OBJ) $(LDFLAGS)

%.o: %.c Makefile $(HEADER)
	@echo [${CC}] $<
	@$(CC) -DSOL_$(SOLSTRUCT) $(CFLAGS) -c $< 

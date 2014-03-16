#######################
#  fltk gui makefile  #
#######################

## build options
# optimized: -Ofast -march=native
# debugging: -g -Og
BUILDOPT = -ggdb3 -Og
# output file
OUTPUT   = ./bin/bouncer

## code folder
SOURCE = src src/math

##############################################################

# compiler + language standard
CC = g++ -std=c++11 $(BUILDOPT)
# compiler flags
CCFLAGS = -c -Wall -Wextra -pedantic `fltk-config --cxxflags` -Iinc
# linker flags
LFLAGS  = -lpthread `fltk-config --ldflags`

##############################################################

# make pipeline
DIRECTORIES = $(SOURCE)
CXXDIRS = $(foreach dir, $(DIRECTORIES), $(dir)/*.cpp)
CXXMODS = $(wildcard $(CXXDIRS))

# compile each .cpp to .o
.cpp.o:
	$(CC) $(CCFLAGS) $< -o $@

# convert .cpp to .o
CXXOBJS = $(CXXMODS:.cpp=.o)
# resource .rc to .o
CCRES   = $(RESOURCES:.rc=.o)

# link all OBJS using CC and link with LFLAGS, then output to OUTPUT
all: $(CXXOBJS)
	$(CC) $(CXXOBJS) $(LFLAGS) -o $(OUTPUT)

# remove each known .o file, and output
clean:
	$(RM) $(CXXOBJS) *~ $(OUTPUT).*


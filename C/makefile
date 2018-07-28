# Makefile for the cacheSim program
CC = g++
CFLAGS = -Wall -g -c
CLINK = $(CC)
OBJS = main.o HopcroftKarp.o lms_structs.o SMG.o Algorithms.o hungarian.o
RM = rm -f
TARGET = decomp

# Creating the  executable
$(TARGET): $(OBJS)
	$(CLINK) -o $(TARGET) $(OBJS)

# Creating the object files
SMG.o: SMG.c SMG.h
lms_structs.o: lms_structs.c lms_structs.h Resources.h
HopcroftKarp.o: HopcroftKarp.c HopcroftKarp.h lms_structs.h Resources.h
hungarian.o: hungarian.c hungarian.h
Algorithms.o: Algorithms.c Algorithms.h lms_structs.h HopcroftKarp.h Resources.h
main.o: main.c Algorithms.h Resources.h
# Cleaning old files before new make
clean:
	$(RM) $(TARGET) *.o *~ "#"* core.*

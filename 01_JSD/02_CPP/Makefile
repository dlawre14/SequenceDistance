############
# Makefile #
############


#########################################################
# Please change values based on your system environment #
#########################################################

# C++ Compiler
CC = g++

# C++ Compile Flag
CFLAGS = -g -Wall -std=gnu++11

# GCC OpenMP Flag
OMPFLAG = -fopenmp

# Static build
# STFLAG = -static


#########################################################
# Please do NOT change below Makefile variables         #
#########################################################
TARGET = JSD

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(OMPFLAG) -o $(TARGET) $(TARGET).cpp

clean:
	$(RM) $(TARGET)

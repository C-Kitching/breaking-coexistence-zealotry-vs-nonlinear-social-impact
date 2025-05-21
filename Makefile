# source /opt/intel/oneapi/setvars.sh
#CC = icpx
CC = g++
#CFLAGS = -std=c++17 -fopenmp -O3 -Wall -g
CFLAGS = -std=c++17 -fopenmp -O3 -Wall -g -I D:/Users/CKitc/eigen-3.4.0/
#CFLAGS = -std=c++17 -fopenmp -O3 -Wall -g -I "C:/Program Files/eigen-3.4.0/"
LFLAGS = -lstdc++fs

TARGET:= main
PARENT_FOLDER := C++
SRC_DIR:= src
BIN_DIR:= bin

ifeq ($(OS),Windows_NT)
    EXE_EXT:=.exe
    PATH_SEP:=\\
    RM= del
else
    EXE_EXT:=
    PATH_SEP:=/
    RM:= rm -f
endif

all: $(PARENT_FOLDER)$(PATH_SEP)$(BIN_DIR)$(PATH_SEP)$(TARGET)$(EXE_EXT)

$(PARENT_FOLDER)$(PATH_SEP)$(BIN_DIR)$(PATH_SEP)$(TARGET)$(EXE_EXT): $(PARENT_FOLDER)$(PATH_SEP)$(SRC_DIR)$(PATH_SEP)$(TARGET).cpp
	$(CC) $(CFLAGS) -o $@ $< $(LFLAGS)

clean:
	$(RM) $(PARENT_FOLDER)$(PATH_SEP)$(BIN_DIR)$(PATH_SEP)$(TARGET)$(EXE_EXT)

run: $(PARENT_FOLDER)$(PATH_SEP)$(BIN_DIR)$(PATH_SEP)$(TARGET)$(EXE_EXT)
	.$(PATH_SEP)$<

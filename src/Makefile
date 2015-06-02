CC=g++
CFLAGS=-c -Wno-c++11-extensions -I../include
LDFLAGS=-O3
SOURCES=io-utils.cc stringprintf.cc
SOURCES_BIN=compute-neighborhood.cc compute-probabilistic-neighborhood.cc compute-randomly-perturbed-neighborhood.cc
OBJECTS=$(SOURCES:.cc=.o)
OBJECTS_BIN=$(SOURCES_BIN:.cc=.o)
EXECUTABLES=compute-neighborhood compute-probabilistic-neighborhood compute-randomly-perturbed-neighborhood

all: $(SOURCES) $(SOURCES_BIN) $(EXECUTABLES)
    
$(EXECUTABLES): $(OBJECTS) $(OBJECTS_BIN)
	$(CC) $(LDFLAGS) $(OBJECTS) $@.o -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLES)
CC=g++
CFLAGS=-c -I../include
LDFLAGS=-O3
SOURCES=test_vectorspace.cc test_neighborgraph.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLES=test_vectorspace test_neighborgraph

all: $(SOURCES) $(EXECUTABLES)
    
$(EXECUTABLES): $(OBJECTS)
	$(CC) $(LDFLAGS) $@.o -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLES)
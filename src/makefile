CC = g++
CFLAGS = -Wall -c -g 
LDFLAGS = -lm
SOURCES = part3.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = exec

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o
	rm exec

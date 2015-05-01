.PHONY : clean

CC = g++
CFLAGS = -std=c++11 -fPIC -O3 -g -lopenblas # C flags
LDFLAGS = -shared -lopenblas  # linking flags

LIB=gropt
SUFFIX=so
SOURCES = $(shell echo src/*.cpp)
HEADERS = $(shell echo incl/*.h)
OBJECTS = $(SOURCES:.cpp=.o)
INC = -Iincl/
TARGET=$(LIB).$(SUFFIX)

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

src/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(TARGET) : $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

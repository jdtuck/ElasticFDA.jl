.PHONY : clean

ERR = $(shell which icpc>/dev/null; echo $$?)
ifeq "$(ERR)" "0"
    CC = icpc
	CFLAGS =  -fPIC -O3 -g # C flags
	LDFLAGS = -shared -lmkl_rt -lpthread # linking flags
else
    CC = g++
	CFLAGS =  -fPIC -std=c++11 -O3 -g -lopenblas # C flags
	LDFLAGS = -shared -lopenblas # linking flags
endif


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

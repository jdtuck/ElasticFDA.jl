.PHONY : clean

ERR = $(shell which icc >/dev/null; echo $$?)
ifeq "$(ERR)" "0"
    CC = icc
else
    CC = gcc
endif

CFLAGS= -fPIC -Wall -Wextra -O3 -g # C flags
LDFLAGS= -shared   # linking flags

LIB=fdasrsf
SUFFIX=so
SOURCES = $(shell echo *.c)
HEADERS = $(shell echo *.h)
OBJECTS=$(SOURCES:.c=.o)

TARGET=$(LIB).$(SUFFIX)

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

install:
	cp $(TARGET) ../

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

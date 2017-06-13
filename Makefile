HOST := $(shell hostname)

ifeq ($(HOST), triolith1)
	CC := mpiicpc # mpi wrapper
	# CC := g++ # This is the main compiler
	# CC := clang --analyze # and comment out the linker last line for sanity
else
	CC := /usr/local/bin/mpic++ # mpi wrapper
	# CC := g++ # This is the main compiler
	# CC := clang --analyze # and comment out the linker last line for sanity
endif

SRCDIR := src
TARGETDIR := bin
TARGET := run

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS := $(SOURCES:%.cpp=%.o)

ifeq ($(HOST), triolith1)
	CFLAGS := -O3 -g3 -ip -xavx -std=c++14 -Wall -pedantic -I Eigen3/ -gxx-name=g++ -cxxlib=/software/apps/gcc/6.2.0/build01
	LIBS := -cxxlib=/software/apps/gcc/6.2.0/build01
else
	CFLAGS := -O3 -g3 -std=c++11 -Wall -pedantic -I Eigen3/
endif

# LIBS := -lsfml-graphics -lsfml-window -lsfml-system

$(TARGETDIR)/$(TARGET): $(OBJECTS)
	@echo "Linking..."
	@mkdir -p $(TARGETDIR)
	@echo "$(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIBS)"; $(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIBS)

%.o: %.$(SRCEXT)
	@echo "$(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<


.PHONY: clean
clean:
	@echo "$(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)"; $(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)

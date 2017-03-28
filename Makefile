CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
TARGETDIR := bin
TARGET := run

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS := $(SOURCES:%.cpp=%.o)

CFLAGS := -O3 -g3 -std=c++11 -Wall -pedantic -I Eigen3/
LIBS := 
# -lsfml-graphics -lsfml-window -lsfml-system

$(TARGETDIR)/$(TARGET): $(OBJECTS)
	@echo "Linking..."
	@mkdir -p $(TARGETDIR)
	@echo "$(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIBS)"; $(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIBS)

%.o: %.$(SRCEXT)
	@echo "$(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<


.PHONY: clean
clean:
	@echo "$(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)"; $(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)












# HOST := $(shell hostname)

# ifeq ($(HOST), triolith1)
# 	CC := icpc # This is the main compiler
# else
# 	CC := g++ # This is the main compiler
# endif

# # CC := clang --analyze # and comment out the linker last line for sanity
# SRCDIR := src
# TARGETDIR := bin
# TARGET := run

# SRCEXT := cpp
# SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
# OBJECTS := $(SOURCES:%.cpp=%.o)

# ifeq ($(CC), icpc)
# 	CFLAGS := -O2 -ip -xavx -g3 -std=c++11 -Wall -pedantic -I Eigen3/ -cxxlib=/software/apps/gcc/4.9.0/build01/
# else
# 	CFLAGS := -O3 -g3 -std=c++11 -Wall -pedantic -I Eigen3/
# endif

# CFLAGS := -O2 -ip -xavx -g3 -std=c++11 -Wall -pedantic -I Eigen3/ -cxxlib=/software/apps/gcc/4.9.0/build01/

# LIBS := -cxxlib=/software/apps/gcc/4.9.0/build01/
# # -lsfml-graphics -lsfml-window -lsfml-system


# $(TARGETDIR)/$(TARGET): $(OBJECTS)
# 	@echo "Linking..."
# 	@mkdir -p $(TARGETDIR)
# 	@echo "$(CC) -gxx-name=g++.orig $^ -o $(TARGETDIR)/$(TARGET) $(LIBS)"; $(CC) -gxx-name=g++.orig $^ -o $(TARGETDIR)/$(TARGET) $(LIBS)

# # @echo "$(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<
# #-gxx-name=g++.orig myprogram.cpp
# %.o: %.$(SRCEXT)
# 	@echo "$(CC) -o $@ -std=c++11 -cxxlib=/software/apps/gcc/4.9.0/build01/ -I Eigen3/ -gxx-name=g++.orig $<"; $(CC) -o $@ -std=c++11 -cxxlib=/software/apps/gcc/4.9.0/build01/ -I Eigen3/ -gxx-name=g++.orig $<


# .PHONY: clean
# clean:
# 	@echo "$(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)"; $(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)
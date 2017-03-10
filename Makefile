CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
TARGETDIR := bin
TARGET := run

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS := $(SOURCES:%.cpp=%.o)

CFLAGS := -O3 -g3 -std=c++11 -Wall -pedantic
LIBS := -lsfml-graphics -lsfml-window -lsfml-system

$(TARGETDIR)/$(TARGET): $(OBJECTS)
	@echo "Linking..."
	@mkdir -p $(TARGETDIR)
	@echo "$(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIB)"; $(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIB)

%.o: %.$(SRCEXT)
	@echo "$(CC) $(CFLAGS) -c -o $@ $<"; $(CC) $(CFLAGS) -c -o $@ $<

clean:
	@echo "$(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)"; $(RM) $(OBJECTS) $(TARGETDIR)/$(TARGET)

.PHONY: clean
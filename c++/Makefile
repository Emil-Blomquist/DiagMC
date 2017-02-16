# PROPAGATOR_DIR = Propagator
# VERTEX_DIR = Vertex
# FEYNMAN_DIAGRAM_DIR = FeynmanDiagram

# MAKE = make

# objects = main.o

# run: output
# 	@echo "----------"
# 	@./output

# output: $(objects)
# 	g++ $(objects) -o output
# 	-lsfml-graphics -lsfml-window -lsfml-system

# # main.o: main.cpp
# # 	g++ -c main.cpp

# clean:
# 	rm *.o output
# 	cd $(PROPAGATOR_DIR) && rm *.o
# 	cd $(VERTEX_DIR) && rm *.o
# 	cd $(FEYNMAN_DIAGRAM_DIR) && rm *.o

CXX=clang++
CXXFLAGS=-g -std=c++11 -Wall -pedantic
BIN=output
LIBS=-lsfml-graphics -lsfml-window -lsfml-system

# SRC=$(wildcard **/*.cpp)
SRC=$(wildcard *.cpp) $(wildcard */*.cpp) $(wildcard */*/*.cpp)
# SRC=main.cpp
OBJ=$(SRC:%.cpp=%.o)

all: $(OBJ)
	@$(CXX) -o $(BIN) $(LIBS) $^
	@echo "-------"
	@./$(BIN)

%.o: %.c
	$(CXX) $@ -c $<

clean:
	@for file in $(OBJ) ; do \
		echo "removing" $$file"..."; \
		rm $$file ; \
	done

	rm $(BIN)
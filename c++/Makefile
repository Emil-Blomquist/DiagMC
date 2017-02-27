# CXX=clang++
# CXXFLAGS=-g -std=c++11 -Wall -pedantic
# BIN=output
# LIBS=-lsfml-graphics -lsfml-window -lsfml-system

# # SRC=$(wildcard **/*.cpp)
# SRC=$(wildcard *.cpp) $(wildcard */*.cpp) $(wildcard */*/*.cpp)
# # SRC=main.cpp
# OBJ=$(SRC:%.cpp=%.o)

# all: $(OBJ)
# 	@$(CXX) -o $(BIN) $(LIBS) $^
# 	@echo "-------"
# 	@./$(BIN)

# %.o: %.c
# 	$(CXX) $@ -c $<

# clean:
# 	@for file in $(OBJ) ; do \
# 		echo "removing" $$file"..."; \
# 		rm $$file ; \
# 	done

# 	rm $(BIN)

CXX=mpic++
CXXFLAGS=-g -std=c++11 -Wall -pedantic
BIN=output
LIBS=-lsfml-graphics -lsfml-window -lsfml-system

# SRC=$(wildcard **/*.cpp)
SRC=$(wildcard *.cpp) $(wildcard */*.cpp) $(wildcard */*/*.cpp)
# SRC=main.cpp
OBJ=$(SRC:%.cpp=%.o)

all: $(OBJ)
	@$(CXX) -o $(BIN) $^

%.o: %.c
	$(CXX) $@ -c $<

clean:
	@for file in $(OBJ) ; do \
		echo "removing" $$file"..."; \
		rm $$file ; \
	done

	rm $(BIN)
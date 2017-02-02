PROPAGATOR_DIR = Propagator
VERTEX_DIR = Vertex

MAKE = make

objects = main.o \
		  $(PROPAGATOR_DIR)/Propagator.o \
		  $(PROPAGATOR_DIR)/Electron.o \
		  $(PROPAGATOR_DIR)/Phonon.o \
		  $(VERTEX_DIR)/Vertex.o

run: output
	@echo "----------"
	@./output

output: $(objects)
	g++ $(objects) -o output

# main.o: main.cpp
# 	g++ -c main.cpp

clean:
	rm *.o output
	cd $(PROPAGATOR_DIR) && rm *.o
	cd $(VERTEX_DIR) && rm *.o
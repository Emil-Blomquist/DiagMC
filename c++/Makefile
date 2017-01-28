PROPAGATOR_DIR = Propagator
MAKE = make

objects = add.o main.o \
		  $(PROPAGATOR_DIR)/Propagator.o

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
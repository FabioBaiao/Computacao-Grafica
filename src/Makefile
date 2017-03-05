COMPILER=g++
GENERATOR_FILES=generator.cpp

build: $(GENERATOR_FILES)
		$(COMPILER) -std=c++0x $(GENERATOR_FILES) -o generator
		@echo 'Use: ./generator <graph_prim> <value1> ... <valueN> <file>'

clean:
	rm generator

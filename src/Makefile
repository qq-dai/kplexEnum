# all: kplexes
# .PHONY : all
#CONFIG = -O3 -g -lpthread

CC = g++ -std=c++14 
CONFIG = -O3 -fopenmp

objects = algorithms.o parallel_plex.o main.o

.PHONY : clean

kplexes: $(objects)
	$(CC) -o kplexes $(objects) $(CONFIG)

%.o:%.cpp
	$(CC) -c $^ $(CONFIG)

clean:
	rm -f *.o kplexes

# all: allplexes
# .PHONY : all
#CONFIG = -O3 -g -lpthread

CC = g++ "-std=c++11" 
CONFIG = -O3 -g

objects = algorithms.o plexenum.o main.o

.PHONY : clean

allplexes: $(objects)
	$(CC) -o allplexes $(objects) $(CONFIG)

%.o:%.cpp
	$(CC) -c $^ $(CONFIG)

clean:
	rm -f *.o allplexes

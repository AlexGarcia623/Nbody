CC = gcc
CFLAGS = -O0 -g

all: main 

main: main.o
	$(CC) $(CFLAGS) -o $@ $^

read_params: read_params.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f main read_params *.o

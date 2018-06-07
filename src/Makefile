# Environment and compilation flags at Exaflop
#   module load gompi
CC=gcc 
# compilation at MacBook
#CC=gcc-7
CFLAGS=-O3 -std=c11 -fopenmp 
LIBS=-lm 

OBJ1=\
	timestep.o \
	derivatives.o \
	source.o \
	utils.o \
	boundary.o \
	map.o

OBJ=main.o $(OBJ1)

TARGET = ModelagemFletcher.exe

all: $(TARGET)

$(TARGET): 	$(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LIBS)


main.o:	main.c $(OBJ1)
	$(CC) -c $(CFLAGS) main.c

timestep.o:	timestep.c timestep.h derivatives.o source.o map.o
	$(CC) -c $(CFLAGS) timestep.c

derivatives.o:	derivatives.c derivatives.h  map.o source.o
	$(CC) -c $(CFLAGS) derivatives.c

boundary.o:	boundary.c boundary.h  map.o
	$(CC) -c $(CFLAGS) boundary.c

source.o:	source.c source.h
	$(CC) -c $(CFLAGS) source.c

utils.o:	utils.c utils.h map.o source.o
	$(CC) -c $(CFLAGS) utils.c

map.o:	map.c map.h
	$(CC) -c $(CFLAGS) map.c

.SUFFIXES	:	.o .c

.c.o:
	$(CC) $(CFLAGS) $*.c

clean:
	rm -f *.o $(TARGET)

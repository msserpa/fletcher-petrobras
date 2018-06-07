# Environment and compilation flags at Exaflop
CC=icc 
CFLAGS=-std=c11 -qopenmp -qopt-report=5 -Wall -Wextra -xHost
LIBS=-lm 

OBJ1=\
	timestep.o \
	derivatives.o \
	source.o \
	utils.o \
	boundary.o \
	map.o

OBJ=main.no.io.o $(OBJ1)
OBJIO=main.io.o $(OBJ1)

TARGET = ModelagemFletcher.no.io.exe
TARGETIO = ModelagemFletcher.io.exe

all: $(TARGET) $(TARGETIO)

$(TARGET): 	$(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LIBS)

$(TARGETIO):      $(OBJIO)
	$(CC) $(CFLAGS) -o $(TARGETIO) $(OBJIO) $(LIBS)

main.no.io.o:	main.no.io.c $(OBJ1)
	$(CC) -c $(CFLAGS) main.no.io.c

main.io.o:   main.io.c $(OBJ1)
	$(CC) -c $(CFLAGS) main.io.c

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
	rm -f *.o *.optrpt $(TARGET) $(TARGETIO)

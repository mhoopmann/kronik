#Set these variables if needed
C = gcc
CC = g++
FLAGS = -O3 -static -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC
INCLUDE = -I.

#Do not touch these variables
KRONIK = CKronik.o


#Make statements
all : Kronik.cpp $(KRONIK)
	$(CC) $(FLAGS) $(INCLUDE) $(KRONIK) Kronik.cpp -o kronik
	ar rcs libkronik.a $(KRONIK)

clean:
	rm *.o kronik


#Kronik objects
CKronik.o : CKronik.cpp
	$(CC) $(FLAGS) $(INCLUDE) CKronik.cpp -c


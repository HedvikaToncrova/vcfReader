CC = g++
CFLAGS = -Wall -std=c++11 -Wno-c++11-extensions

all: main.o vcfData.o vcfParser.o
	$(CC) main.o vcfData.o vcfParser.o -o vcfReader

test: unitTest.cpp vcfData.o vcfParser.o
	$(CC) $(CFLAGS) unitTest.cpp vcfData.o vcfParser.o -lboost_unit_test_framework -o unitTest 


main.o: main.cpp vcfData.hpp
	$(CC) $(CFLAGS) -c main.cpp

vcfData.o: vcfData.cpp vcfData.hpp vcfParser.hpp
	$(CC) $(CFLAGS) -c vcfData.cpp

vcfParser.o: vcfParser.cpp vcfParser.hpp
	$(CC) $(CFLAGS) -c vcfParser.cpp

clean:
	rm -rf *o vcfReader unitTest

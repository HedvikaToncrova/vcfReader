CC = g++
CFLAGS = -Wall -c -std=c++11 -Wno-c++11-extensions

all: main.o vcfData.o vcfParser.o
	$(CC) main.o vcfData.o vcfParser.o -o vcfReader

main.o: main.cpp vcfData.o
	$(CC) $(CFLAGS) main.cpp

vcfData.o: vcfData.cpp vcfData.hpp vcfParser.hpp
	$(CC) $(CFLAGS) vcfData.cpp

vcfParser.o: vcfParser.cpp vcfParser.hpp
	$(CC) $(CFLAGS) vcfParser.cpp



test: unitTest.cpp vcfData.hpp
	$(CC) $(CFLAGS)  unitTest.cpp vcfData.cpp -lboost_unit_test_framework -o unitTest

clean:
	rm -rf *o vcfReader unitTest

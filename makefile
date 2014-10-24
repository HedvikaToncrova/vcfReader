CC = g++
CFLAGS = -Wall -std=c++11 -Wno-c++11-extensions

all: main.cpp vcfData.hpp
	$(CC) $(CFLAGS) main.cpp vcfData.cpp -o vcfReader

test: unitTest.cpp vcfData.hpp
	$(CC) $(CFLAGS)  unitTest.cpp vcfData.cpp -lboost_unit_test_framework -o unitTest

clean:
	rm -rf *o vcfReader unitTest

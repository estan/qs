CXXFLAGS=-O2 -std=gnu++0x -Wall -pedantic
LDFLAGS=-O2
LIBS=-lgmpxx -lgmp

all: factor

factor: main.cpp fermat.o rho.o qs.o
	g++ -c $(CXXFLAGS) -o main.o main.cpp
	g++ $(LDFLAGS) -o factor main.o fermat.o rho.o qs.o $(LIBS)

rho.o: rho.cpp rho.h
	g++ -c $(CXXFLAGS) -o rho.o rho.cpp

fermat.o: fermat.cpp fermat.h
	g++ -c $(CXXFLAGS) -o fermat.o fermat.cpp

qs.o: qs.cpp qs.h matrix.h
	g++ -c $(CXXFLAGS) -o qs.o qs.cpp

clean:
	-@rm *.o factor 2> /dev/null || true

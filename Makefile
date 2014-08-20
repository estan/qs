CXXFLAGS=-O2 -std=gnu++0x -Wall -pedantic
LDFLAGS=-O2
LIBS=-lgmpxx -lgmp

all: factor

factor: main.cpp fermat.o rho.o qs.o
	$(CXX) -c $(CXXFLAGS) -o main.o main.cpp
	$(CXX) $(LDFLAGS) -o factor main.o fermat.o rho.o qs.o $(LIBS)

rho.o: rho.cpp rho.h
	$(CXX) -c $(CXXFLAGS) -o rho.o rho.cpp

fermat.o: fermat.cpp fermat.h
	$(CXX) -c $(CXXFLAGS) -o fermat.o fermat.cpp

qs.o: qs.cpp qs.h matrix.h
	$(CXX) -c $(CXXFLAGS) -o qs.o qs.cpp

clean:
	-@rm *.o factor 2> /dev/null || true

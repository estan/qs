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

test: test-50bits test-55bits test-60bits test-65bits test-70bits test-75bits \
	test-80bits test-85bits test-90bits test-95bits test-100bits

test-%: factor
	@echo running tests in data/test$*.in
	@while read N p q; do \
		actual=$$(echo $$N | ./factor | sort); \
		expected=$$(echo -e "$$p\n$$q\n" | sort); \
		if [[ $$actual != $$expected ]]; then \
		   echo test case $$N $$p $$q in data/test$*.in failed; \
		   echo expected: $$expected, but got $$actual; \
		   exit 1; \
		fi \
	done < data/test$*.in

clean:
	-@rm *.o factor 2> /dev/null || true

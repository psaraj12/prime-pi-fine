CXX = g++
CXXFLAGS = -O3 -march=native -fopenmp -std=c++17

TARGET = prime_pi_fine
SRC = src/prime_pi_fine_autoalign.cpp

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
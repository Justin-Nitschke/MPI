CXX = mpicxx

CXXFLAGS = -Wall -Wextra -O3

TARGET = main

SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	mpiexec -n 4 ./$(TARGET) -i in/test.txt

clean:
	rm -f $(TARGET)
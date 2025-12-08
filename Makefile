CXX = mpicxx

CXXFLAGS = -Wall -Wextra -O3

LDFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw3 -lX11 -lXrandr -lpthread -ldl

TARGET = nbody

SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	mpiexec -n 4 ./$(TARGET) -i in/test.txt -o out/results.txt -t  0.5 -s 2000 -d 0.005 -v

clean:
	rm -f $(TARGET)
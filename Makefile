CXX = mpicxx

CXXFLAGS = -Wall -Wextra -O3

LDFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -ldl

TARGET = main

SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	mpiexec -n 4 ./$(TARGET) -i in/test.txt -o out/result.txt -v

clean:
	rm -f $(TARGET)
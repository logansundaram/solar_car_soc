CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra
LDFLAGS :=

SRC := main.cpp Battery.cpp Battery_Test.cpp Pack_test.cpp
OBJ := $(SRC:.cpp=.o)

TARGET := main.exe

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDFLAGS)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


release: CXXFLAGS := -std=c++17 -O3 -DNDEBUG
release: clean $(TARGET)
	@echo "Built release version."

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all release clean

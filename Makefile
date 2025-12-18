CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native



TARGET2 = test_fast_put
SOURCE2 = test_fast_put.cpp
TARGET3 = reference_tree_put
SOURCE3 = reference_tree_put.cpp

all:  $(TARGET2) $(TARGET3)


$(TARGET2): $(SOURCE2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(SOURCE2)

$(TARGET3): $(SOURCE3)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(SOURCE3)

clean:
	rm -f  $(TARGET2) $(TARGET3)

.PHONY: clean all


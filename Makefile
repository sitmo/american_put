CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native



TARGET2 = americanput
SOURCE2 = americanput.cpp
TARGET3 = americanput_tree
SOURCE3 = americanput_tree.cpp

all:  $(TARGET2) $(TARGET3)


$(TARGET2): $(SOURCE2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(SOURCE2)

$(TARGET3): $(SOURCE3)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(SOURCE3)

clean:
	rm -f  $(TARGET2) $(TARGET3)

.PHONY: clean all


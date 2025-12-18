CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native



TARGET2 = test_fast_put
SOURCE2 = test_fast_put.cpp
TARGET3 = reference_tree_put
SOURCE3 = reference_tree_put.cpp
TARGET4 = test_table2
SOURCE4 = test_table2.cpp
TARGET5 = test_table3
SOURCE5 = test_table3.cpp
TARGET6 = test_ref33_table2
SOURCE6 = test_ref33_table2.cpp
TARGET7 = test_ref33_table3
SOURCE7 = test_ref33_table3.cpp

all:  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)


$(TARGET2): $(SOURCE2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(SOURCE2)

$(TARGET3): $(SOURCE3)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(SOURCE3)

$(TARGET4): $(SOURCE4)
	$(CXX) $(CXXFLAGS) -o $(TARGET4) $(SOURCE4)

$(TARGET5): $(SOURCE5)
	$(CXX) $(CXXFLAGS) -o $(TARGET5) $(SOURCE5)

$(TARGET6): $(SOURCE6)
	$(CXX) $(CXXFLAGS) -o $(TARGET6) $(SOURCE6)

$(TARGET7): $(SOURCE7)
	$(CXX) $(CXXFLAGS) -o $(TARGET7) $(SOURCE7)

clean:
	rm -f  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)

.PHONY: clean all


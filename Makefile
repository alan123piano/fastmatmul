EXECUTABLE 	= fastmatmul
CXX         = g++
CXXFLAGS 	= -std=c++1z -Wconversion -Wall -Werror -Wextra -pedantic -O3 -DNDEBUG

all: $(EXECUTABLE)

$(EXECUTABLE): $(EXECUTABLE).o
	$(CXX) $(CXXFLAGS) $(EXECUTABLE).o -o $(EXECUTABLE)
	rm $(EXECUTABLE).o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp

clean:
	rm -f $(EXECUTABLE) *.o


CXX = g++
CXXFLAGS = -Wall -g

HTSLIB_DIR = ../htslib

FERMILIB_DIR = ../fermi-lite

EDLIB_DIR = ../edlib/include

INCLUDES = -I $(HTSLIB_DIR) -I $(FERMILIB_DIR) -I $(EDLIB_DIR) -fopenmp
LIBS = -L$(HTSLIB_DIR) -L$(FERMILIB_DIR) -lfml -lz -lbz2 -llzma -lpthread -lhts -fopenmp

# Source and object files
SRCS = rivindel.cpp BamReader.cpp BamRecord.cpp BamWriter.cpp RefFasta.cpp VcfWriter.cpp ../edlib/src/edlib.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

TARGET = rivindel

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

-include $(DEPS)

.PHONY: clean run

clean:
	rm -f $(OBJS) $(TARGET) $(DEPS)

run: $(TARGET)
	./$(TARGET)



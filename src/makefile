# Assume these variables are defined earlier in your Makefile
CXX = g++
CXXFLAGS = -Wall -g -std=c++17
HTSLIB_DIR = ../htslib
FERMILIB_DIR = ../fermi-lite
EDLIB_DIR = ../edlib/include
SSW_DIR = .

# Include directories and libraries
INCLUDES = -I$(HTSLIB_DIR) -I$(FERMILIB_DIR) -I$(EDLIB_DIR) -I$(SSW_DIR) -fopenmp
LIBS = -L$(HTSLIB_DIR) -L$(FERMILIB_DIR) -lfml -lz -lbz2 -llzma -lpthread -lhts -fopenmp

# Source files
SRCS = rivindel.cpp BamReader.cpp BamRecord.cpp BamWriter.cpp RefFasta.cpp VcfReader.cpp VcfWriter.cpp Assembler.cpp ../edlib/src/edlib.cpp $(SSW_DIR)/ssw.cpp $(SSW_DIR)/ssw_cpp.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
OBJS += $(SSW_DIR)/ssw_cpp.o

# Target executable
TARGET = rivindel

# Compilation rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


# Linking rule for the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Include dependency files
-include $(DEPS)

# Phony targets
.PHONY: all clean run

# Default target
all: $(TARGET)

# Clean target
clean:
	$(RM) -f $(OBJS) $(TARGET) $(DEPS)

# Run target
run: $(TARGET)
	./$(TARGET)



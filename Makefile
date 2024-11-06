SHELL = bash

COMPILER = g++
CXXFLAGS = -O3 -shared -std=c++17 -fPIC -flto -Wall -Wextra -Wattributes

PROJECT=depmod

PY = "3.11"
PYTHON_INC = $(shell python$(PY)-config --includes)
PYTHON_EXT = $(shell python$(PY)-config --extension-suffix)

EIGEN = external/eigen
PYBIND11 = $(shell python$(PY) -c "import pybind11; print(pybind11.get_include())")

PROJECTDIR = $(realpath $(CURDIR))
SOURCESDIR = $(PROJECTDIR)/depmod-lib


INCLUDES = -I$(SOURCESDIR)/include $(PYTHON_INC) -I$(PYBIND11) -I$(EIGEN)

SRC  = $(wildcard $(SOURCESDIR)/src/*.cpp) $(wildcard $(SOURCESDIR)/binding/*.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = $(SRC:%.cpp=%.d)


LIB=_lib$(PYTHON_EXT)

all: $(LIB)
	$(info ${SRC})

$(LIB): $(OBJS)
	$(COMPILER) $(CXXFLAGS) $(INCLUDES) -o $@ $+ -shared

-include $(DEPS)

%.o: %.cpp Makefile
	$(COMPILER) $(CXXFLAGS) $(INCLUDES) -MMD -MP -o $@ -c $<

clean:
	rm -rf $(LIB)
	rm -rf $(OBJS)
	rm -rf $(DEPS)
	rm -rf $(PROJECT)/$(LIB)

install:
	mv $(LIB) $(PROJECT)/$(LIB)

test:
	$(SHELL) $(PROJECTDIR)/scripts/run_tests.sh

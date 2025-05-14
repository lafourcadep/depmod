SHELL = bash

PROJECT=depmod

COMPILER = g++
CXXFLAGS := -O3 -shared -std=c++20 -fPIC -flto=auto -Wall -Wextra -Wattributes
LDFLAGS :=

PY = "3.11"
PYTHON_INC = $(shell python$(PY)-config --includes)
PYTHON_EXT = $(shell python$(PY)-config --extension-suffix)
EIGEN = external/eigen

ZLIB_OK := $(shell ./scripts/check_zlib.sh >/dev/null 2>&1 && echo yes || echo no)
BZIP2_OK := $(shell ./scripts/check_bzip2.sh >/dev/null 2>&1 && echo yes || echo no)
LZMA_OK := $(shell ./scripts/check_lzma.sh >/dev/null 2>&1 && echo yes || echo no)

ifeq ($(ZLIB_OK),yes)
  CXXFLAGS += -DUSE_ZLIB
  LDFLAGS += -lz
endif

ifeq ($(BZIP2_OK),yes)
  CXXFLAGS += -DUSE_BZIP2
  LDFLAGS += -lbz2
endif

ifeq ($(LZMA_OK),yes)
  CXXFLAGS += -DUSE_LZMA
  LDFLAGS += -llzma
endif

PYBIND11 = $(shell python$(PY) -c "import pybind11; print(pybind11.get_include())")
PROJECTDIR = $(realpath $(CURDIR))
SOURCESDIR = $(PROJECTDIR)/src/cpp

INCLUDES = -I$(SOURCESDIR) $(PYTHON_INC) -I$(PYBIND11) -I$(EIGEN)

SRC  = $(wildcard $(SOURCESDIR)/*.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = $(SRC:%.cpp=%.d)

LIB=_lib$(PYTHON_EXT)

all: $(LIB)
	$(info ${SRC})

$(LIB): $(OBJS)
	$(COMPILER) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $+ -shared 

-include $(DEPS)

%.o: %.cpp Makefile
	$(COMPILER) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -MMD -MP -o $@ -c $<

clean:
	rm -rf $(LIB)
	rm -rf $(OBJS)
	rm -rf $(DEPS)
	rm -rf $(PROJECTDIR)/src/$(PROJECT)/$(LIB)

install:
	mv $(LIB) $(PROJECTDIR)/src/$(PROJECT)/$(LIB)

test:
	$(SHELL) $(PROJECTDIR)/scripts/run_tests.sh

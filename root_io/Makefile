ROOTCONFIG ?= root-config

CXX      := $(shell $(ROOTCONFIG) --cxx)
CXXFLAGS := -O2 -Wall -Wextra $(shell $(ROOTCONFIG) --cflags)
LDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
LIBS     := $(shell $(ROOTCONFIG) --libs)

all: clumpy2root

clumpy2root: clumpy2root.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -f clumpy2root

.PHONY: all clean

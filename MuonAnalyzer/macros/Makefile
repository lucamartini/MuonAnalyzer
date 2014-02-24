ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -O2 -Wall -Wextra 
LDFLAGS       = -g3

SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

CXXFLAGS      += $(ROOTCFLAGS)
NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../../../../../lib/

###########
# TARGETS #
###########

mass: mass.cc 
		$(CXX) $(CXXFLAGS) -o mass $(GLIBS)  $(LDFLAGS)  $<


clean:
	rm mass



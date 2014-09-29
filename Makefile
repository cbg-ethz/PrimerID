export CXX = clang++
export CXXFLAGS = -O3 -Wall -pedantic -std=gnu++11
export CPPFLAGS = -I. -I/opt/local/include/
export LDFLAGS = -L/opt/local/lib -Wl,-dead_strip_dylibs
export LIBS = -lgsl

SUBDIRS = src

.PHONY: clean default clean cleanobj

default: PID

clean: cleanobj
	rm -rf PID

cleanobj:
	rm -rf src/*.o

PID: subdirs

subdirs:
	$(MAKE) -C $(SUBDIRS)

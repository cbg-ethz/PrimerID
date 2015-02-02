ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = pidalyser

pidalyser_includes = \
	src/alignment.hpp \
	src/kmeans.hpp \
	src/prob_cycle.hpp \
	src/proper_read.hpp \
	src/reference.hpp \
	src/statistics.hpp
pidalyser_SOURCES = $(pidalyser_includes) \
	src/alignment.cpp \
	src/main.cpp \
	src/proper_read.cpp \
	src/reference.cpp
pidalyser_CFLAGS = $(GSL_CFLAGS)
pidalyser_LDADD = $(GSL_LIBS)
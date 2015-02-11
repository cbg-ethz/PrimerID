ACLOCAL_AMFLAGS = -I m4

MAINTAINERCLEANFILES = Makefile Makefile.in \
	aclocal.m4 autom4te.cache/ \
	config.h config.h.in config.log config.status configure \
	depcomp install-sh missing \
	pidalyser src/.deps/ src/.libs/ src/.dirstamp stamp-h1

bin_PROGRAMS = pidalyser pidalign

# pidalyser
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

# pidalign
pidalign_includes = \
	src/ednafull.hpp
pidalign_SOURCES = $(pidalign_includes) \
	src/pidalign.cpp
pidalign_CPPFLAGS = \
	-I$(SEQAN_INCLUDEDIR)
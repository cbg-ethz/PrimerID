ACLOCAL_AMFLAGS = -I m4

MAINTAINERCLEANFILES = Makefile Makefile.in \
	aclocal.m4 autom4te.cache/ \
	config.h config.h.in config.log config.status configure \
	depcomp install-sh missing \
	pidalyse src/.deps/ src/.libs/ src/.dirstamp stamp-h1

bin_PROGRAMS = pidalyse pidalign

# pidalyse
pidalyse_includes = \
	src/alignment.hpp \
	src/kmeans.hpp \
	src/prob_cycle.hpp \
	src/proper_read.hpp \
	src/reference.hpp \
	src/statistics.hpp
pidalyse_SOURCES = $(pidalyse_includes) \
	src/alignment.cpp \
	src/main.cpp \
	src/proper_read.cpp \
	src/reference.cpp
pidalyse_CXXFLAGS = $(GSL_CFLAGS)
pidalyse_LDADD = $(GSL_LIBS)

# pidalign
libthreadpool11_sources = \
	src/threadpool11/pool.hpp \
	src/threadpool11/threadpool11.hpp \
	src/threadpool11/utility.hpp \
	src/threadpool11/worker.hpp \
	src/threadpool11/pool.cpp \
	src/threadpool11/worker.cpp

pidalign_includes = \
	src/ednafull.hpp
pidalign_SOURCES = $(pidalign_includes) $(libthreadpool11_sources) \
	src/pidalign.cpp
pidalign_CPPFLAGS = \
	-I$(SEQAN_INCLUDEDIR) $(BOOST_CPPFLAGS)
pidalign_LDFLAGS = \
	$(BOOST_LDFLAGS)
bin_PROGRAMS = pidalyse pidalign

# pidalyse
pidalyse_includes = \
	src/pidalyse/alignment.hpp \
	src/pidalyse/kmeans.hpp \
	src/pidalyse/prob_cycle.hpp \
	src/pidalyse/proper_read.hpp \
	src/pidalyse/reference.hpp \
	src/pidalyse/statistics.hpp
pidalyse_SOURCES = $(pidalyse_includes) \
	src/pidalyse/alignment.cpp \
	src/pidalyse/main.cpp \
	src/pidalyse/proper_read.cpp \
	src/pidalyse/reference.cpp
pidalyse_CXXFLAGS = $(GSL_CFLAGS)
pidalyse_LDADD = $(GSL_LIBS)

# pidalign
libthreadpool11_sources = \
	src/pidalign/threadpool11/pool.hpp \
	src/pidalign/threadpool11/threadpool11.hpp \
	src/pidalign/threadpool11/utility.hpp \
	src/pidalign/threadpool11/worker.hpp \
	src/pidalign/threadpool11/pool.cpp \
	src/pidalign/threadpool11/worker.cpp

pidalign_includes = \
	src/pidalign/ednafull.hpp
pidalign_SOURCES = $(pidalign_includes) $(libthreadpool11_sources) \
	src/pidalign/pidalign.cpp
pidalign_CPPFLAGS = \
	-I$(SEQAN_INCLUDEDIR) $(BOOST_CPPFLAGS)
pidalign_LDFLAGS = \
	$(BOOST_LDFLAGS)

bin_PROGRAMS = pidalyse pidalign pidrtpcrsim

# pidalyse
pidalyse_includes = \
	src/pidalyse/alignment.hpp \
	src/pidalyse/proper_read.hpp \
	src/pidalyse/reference.hpp \
	src/pidalyse/statistics.hpp
pidalyse_SOURCES = $(pidalyse_includes) \
	src/pidalyse/alignment.cpp \
	src/pidalyse/main.cpp \
	src/pidalyse/proper_read.cpp \
	src/pidalyse/reference.cpp \
	src/pidalyse/statistics.cpp
pidalyse_CXXFLAGS = $(GSL_CFLAGS)
pidalyse_LDADD = $(GSL_LIBS)  $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)


# pidalign
boost_threadpool_sources = \
	src/pidalign/threadpool/future.hpp \
	src/pidalign/threadpool/pool.hpp \
	src/pidalign/threadpool/pool_adaptors.hpp \
	src/pidalign/threadpool/scheduling_policies.hpp \
	src/pidalign/threadpool/shutdown_policies.hpp \
	src/pidalign/threadpool/size_policies.hpp \
	src/pidalign/threadpool/task_adaptors.hpp \
	src/pidalign/threadpool/threadpool.hpp \
	src/pidalign/threadpool/detail/future.hpp \
	src/pidalign/threadpool/detail/locking_ptr.hpp \
	src/pidalign/threadpool/detail/pool_core.hpp \
	src/pidalign/threadpool/detail/scope_guard.hpp \
	src/pidalign/threadpool/detail/worker_thread.hpp
pidalign_includes = \
	src/pidalign/ednafull.hpp
pidalign_SOURCES = $(pidalign_includes) $(boost_threadpool_sources) \
	src/pidalign/pidalign.cpp
pidalign_CPPFLAGS = \
	-I$(SEQAN_INCLUDEDIR) $(BOOST_CPPFLAGS)
pidalign_LDFLAGS = \
	$(BOOST_LDFLAGS)
pidalign_LDADD = $(BOOST_THREAD_LIB) $(BOOST_SYSTEM_LIB)


# pidrtpcrsim
pidrtpcrsim_SOURCES = src/pidrtpcrsim/pidrtpcrsim.cpp
pidrtpcrsim_CPPFLAGS = $(GSL_CFLAGS) $(BOOST_CPPFLAGS)
pidrtpcrsim_LDFLAGS = $(BOOST_LDFLAGS)
pidrtpcrsim_LDADD = $(GSL_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB) $(BOOST_SYSTEM_LIB)



# CHECKS
check_PROGRAMS = test-dnavector test-prob-cycle
test_dnavector_SOURCES = test/test_dnavector.cpp src/pidalyse/statistics.cpp src/pidalyse/proper_read.cpp
test_dnavector_CPPFLAGS = $(BOOST_CPPFLAGS)
test_dnavector_LDFLAGS = $(GSL_LIBS)

test_prob_cycle_SOURCES = test/test_prob_cycle.cpp src/pidalyse/statistics.cpp src/pidalyse/proper_read.cpp
test_prob_cycle_CPPFLAGS = $(BOOST_CPPFLAGS)
test_prob_cycle_LDFLAGS = $(GSL_LIBS)
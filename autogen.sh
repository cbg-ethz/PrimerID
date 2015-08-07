#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	# executables
	rm -rf pidalign
	rm -rf pidalyse
	rm -rf pidrtpcrsim
	rm -rf test-dnavector
	rm -rf test-prob-cycle

	# build associated dependencies
	rm -rf src/pidalign/.deps/
	rm -rf src/pidalign/.dirstamp
	rm -rf src/pidalign/pidalign-pidalign.o
	
	rm -rf src/pidalyse/.deps/
	rm -rf src/pidalyse/.dirstamp
	rm -rf src/pidalyse/pidalyse-alignment.o
	rm -rf src/pidalyse/pidalyse-main.o
	rm -rf src/pidalyse/pidalyse-proper_read.o
	rm -rf src/pidalyse/pidalyse-reference.o
	rm -rf src/pidalyse/pidalyse-statistics.o
	
	rm -rf src/pidrtpcrsim/.deps/
	rm -rf src/pidrtpcrsim/.dirstamp
	rm -rf src/pidrtpcrsim/pidrtpcrsim-pidrtpcrsim.o
	
	rm -rf src/pidalyse/test_dnavector-proper_read.o
	rm -rf src/pidalyse/test_dnavector-statistics.o
	rm -rf src/pidalyse/test_prob_cycle-proper_read.o
	rm -rf src/pidalyse/test_prob_cycle-statistics.o
	
	rm -rf test/.deps/
	rm -rf test/.dirstamp
	rm -rf test/test_dnavector-test_dnavector.o
	rm -rf test/test_prob_cycle-test_prob_cycle.o
	
	# autotools cruft
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf compile
	rm -rf config.guess
	rm -rf config.h
	rm -rf config.h.in
	rm -rf config.h.in~
	rm -rf config.log
	rm -rf config.status
	rm -rf config.sub
	rm -rf configure
	rm -rf depcomp
	rm -rf install-sh
	rm -rf missing
	rm -rf stamp-h1
	
	# OS X cruft
	find . -name '.DS_Store' -type f -delete
	
	# R cruft
	find . -name '.RData' -type f -delete
	find . -name '.Rapp.history' -type f -delete
}

if [[ "$1" == "--clean" ]]
then
	clean_files
	exit
fi

echo "Bootstrapping Autotools"
autoreconf -vif

if [[ "$1" == "--test" ]]
then
	echo "${DISTCHECK_CONFIGURE_FLAGS}"
	./configure ${DISTCHECK_CONFIGURE_FLAGS}
	DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS}" make distcheck
	clean_files
fi

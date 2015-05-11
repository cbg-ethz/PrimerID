#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	# ./
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

	# ./src/pidalign
	rm -rf pidalign
	rm -rf src/pidalign/.deps/
	rm -rf src/pidalign/.dirstamp
	rm -rf src/pidalign/*.o

	# ./src/pidalyse
	rm -rf pidalyse
	rm -rf src/pidalyse/.deps/
	rm -rf src/pidalyse/.dirstamp
	rm -rf src/pidalyse/*.o
	
	# ./src/pidrtpcrsim
	rm -rf pidrtpcrsim
	rm -rf src/pidrtpcrsim/.deps/
	rm -rf src/pidrtpcrsim/.dirstamp
	rm -rf src/pidrtpcrsim/*.o
	
	# ./test
	rm -rf test-dnavector
	rm -rf test-prob-cycle
	rm -rf test/.deps/
	rm -rf test/.dirstamp
	rm -rf test/*.o

	# Tarballs
	rm -rf pidalyse-0.1.tar.bz2
	rm -rf pidalyse-0.1/

	# OS X files
	rm -rf .DS_Store
	rm -rf src/.DS_Store
	rm -rf src/pidalign/.DS_Store
	rm -rf src/pidalyse/.DS_Store
	rm -rf src/pidrtpcrsim/.DS_Store
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

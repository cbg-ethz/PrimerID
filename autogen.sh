#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	# ./
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf config.h
	rm -rf config.h.in
	rm -rf config.h.in~
	rm -rf config.log
	rm -rf config.guess
	rm -rf config.status
	rm -rf config.sub
	rm -rf configure
	rm -rf depcomp
	rm -rf install-sh
	rm -rf missing
	rm -rf stamp-h1

	rm -rf pidalign
	rm -rf pidalyse

	# ./src
	rm -rf src/.deps/
	rm -rf src/.dirstamp
	rm -rf src/.libs/
	rm -rf src/alignment.o
	rm -rf src/main.o
	rm -rf src/pidalign-pidalign.o
	rm -rf src/proper_read.o
	rm -rf src/reference.o

	rm -rf src/threadpool11/.deps/
	rm -rf src/threadpool11/.dirstamp
	rm -rf src/threadpool11/.libs/
	rm -rf src/threadpool11/pool.o
	rm -rf src/threadpool11/worker.o

	# Tarballs
	rm -rf pidalyse-0.1.tar.bz2
	rm -rf pidalyse-0.1/

	# OS X files
	rm -rf .DS_Store src/.DS_Store m4/.DS_Store src/threadpool11/.DS_Store
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

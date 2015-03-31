#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	# ./
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf config.guess
	rm -rf config.h
	rm -rf config.h.in
	rm -rf config.log
	rm -rf config.status
	rm -rf config.sub
	rm -rf configure
	rm -rf depcomp
	rm -rf install-sh
	rm -rf missing
	rm -rf pidalign
	rm -rf pidalyse
	rm -rf stamp-h1

	# ./src/pidalign
	rm -rf src/pidalign/.deps/
	rm -rf src/pidalign/.dirstamp
	rm -rf src/pidalign/pidalign-pidalign.o
	rm -rf src/pidalign/threadpool11/.deps/
	rm -rf src/pidalign/threadpool11/.dirstamp
	rm -rf src/pidalign/threadpool11/pidalign-pool.o
	rm -rf src/pidalign/threadpool11/pidalign-worker.o

	# ./src/pidalyse
	rm -rf src/pidalyse/.deps/
	rm -rf src/pidalyse/.dirstamp
	rm -rf src/pidalyse/pidalyse-alignment.o
	rm -rf src/pidalyse/pidalyse-main.o
	rm -rf src/pidalyse/pidalyse-proper_read.o
	rm -rf src/pidalyse/pidalyse-reference.o

	# Tarballs
	rm -rf pidalyse-0.1.tar.bz2
	rm -rf pidalyse-0.1/

	# OS X files
	rm -rf .DS_Store
	rm -rf src/.DS_Store
	rm -rf m4/.DS_Store
	rm -rf src/pidalign/.DS_Store
	rm -rf src/pidalign/threadpool11/.DS_Store
	rm -rf src/pidalyse/.DS_Store
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

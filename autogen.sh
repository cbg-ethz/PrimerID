#!/bin/sh

if [[ "$1" == "--clean" ]]
then
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
	rm -rf config.status
	rm -rf configure
	rm -rf depcomp
	rm -rf install-sh
	rm -rf missing
	rm -rf stamp-h1

	rm -rf pidalign
	rm -rf pidalyser

	# ./src
	rm -rf src/.deps/
	rm -rf src/.dirstamp
	rm -rf src/.libs/
	rm -rf src/alignment.o
	rm -rf src/main.o
	rm -rf src/pidalign-pidalign.o
	rm -rf src/proper_read.o
	rm -rf src/reference.o

	# Tarballs
	rm -rf pidalyser-0.1.tar.bz2
	rm -rf pidalyser-0.1/

	# OS X files
	rm -rf .DS_Store src/.DS_Store m4/.DS_Store
else
	echo "Bootstrapping Autotools"
	autoreconf -vif
fi

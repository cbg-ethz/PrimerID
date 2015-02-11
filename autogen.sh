#!/bin/sh

if [[ "$1" == "--clean" ]]
then
	echo "Cleaning bootstrapped files"

	# ./
	rm -rf tmp/
	rm -rf .deps/
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf autoscan-2.69.log
	rm -rf compile
	rm -rf config.guess
	rm -rf config.h
	rm -rf config.h.in
	rm -rf config.h.in~
	rm -rf config.log
	rm -rf config.status
	rm -rf config.sub
	rm -rf configure
	rm -rf configure.scan
	rm -rf depcomp
	rm -rf install-sh
	rm -rf libtool
	rm -rf ltmain.sh
	rm -rf m4/libtool.m4
	rm -rf m4/ltoptions.m4
	rm -rf m4/ltsugar.m4
	rm -rf m4/ltversion.m4
	rm -rf m4/lt~obsolete.m4
	rm -rf missing
	rm -rf stamp-h1
	rm -rf LogLik.R
	rm -rf LogLik.pdf
	rm -rf LogLik.Rout
	rm -rf .RData
	rm -rf .Rapp.history
	rm -rf pidalyser
	rm -rf pidalign

	# ./src
	rm -rf src/.dirstamp
	rm -rf src/.deps/
	rm -rf src/Makefile
	rm -rf src/Makefile.in
	rm -rf src/alignment.o
	rm -rf src/main.o
	rm -rf src/proper_read.o
	rm -rf src/reference.o
	rm -rf src/pidalign-pidalign.o

	# Tarballs
	rm -rf pidalyser-0.1.tar.bz2
	rm -rf pidalyser-0.1.tar.gz
	rm -rf pidalyser-0.1/

	# OS X files
	rm -rf .DS_Store src/.DS_Store
else
	echo "Bootstrapping Autotools"
	autoreconf -vif
fi

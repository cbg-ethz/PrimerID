dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[pIDalyse],
	[0.1],
	[david.seifert@bsse.ethz.ch],
	[pidalyse],
	[https://github.com/cbg-ethz/PrimerID])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/pidalyse/main.cpp])
AC_LANG([C++])



dnl ========================================================================
dnl Check whether we want to set defaults for CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ========================================================================
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [Modifies various CXXFLAGS, CPPFLAGS and LDFLAGS required for debugging])]
)

AS_IF(
	[test "x$enable_debug" = "xyes"],
		[
			enable_cxxflags_setting=yes
			cxxflags_test="-O0 -g"
			
			enable_cppflags_setting=yes
			cppflags_test=""
			
			enable_ldflags_setting=yes
			ldflags_test=""
		],
		[
			cxxflags_test="-O2"
			cppflags_test="-DNDEBUG -DSEQAN_ENABLE_TESTING=0"
			ldflags_test="-Wl,-O1 -Wl,--as-needed -Wl,-dead_strip_dylibs"
		
			AC_MSG_CHECKING([whether configure should try to set CXXFLAGS])
			AS_IF(
				[test "x${CXXFLAGS+set}" = "xset"],
					[enable_cxxflags_setting=no],
					[enable_cxxflags_setting=yes]
			)
			AC_MSG_RESULT([$enable_cxxflags_setting])
                	
			AC_MSG_CHECKING([whether configure should try to set CPPFLAGS])
			AS_IF(
				[test "x${CPPFLAGS+set}" = "xset"],
					[enable_cppflags_setting=no],
					[enable_cppflags_setting=yes]
			)
			AC_MSG_RESULT([$enable_cppflags_setting])
                	
			AC_MSG_CHECKING([whether configure should try to set LDFLAGS])
			AS_IF(
				[test "x${LDFLAGS+set}" = "xset"],
					[enable_ldflags_setting=no],
					[enable_ldflags_setting=yes]
			)
			AC_MSG_RESULT([$enable_ldflags_setting])
		]
)



dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AX_COMPILER_VENDOR
AX_COMPILER_VERSION
AC_PROG_SED
AC_PROG_CXX
AC_CHECK_HEADERS([cstdint])



dnl ==================================
dnl Set CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ==================================
AS_IF(
	[test "x${enable_cxxflags_setting}" = "xyes"],
		[
		CXXFLAGS=""
		AX_APPEND_COMPILE_FLAGS([$cxxflags_test -Wall -Wno-redeclared-class-member], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_cppflags_setting}" = "xyes"],
		[
		CPPFLAGS=""
		AX_APPEND_COMPILE_FLAGS([$cppflags_test], [CPPFLAGS])
		CPPFLAGS=$( echo $CPPFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_ldflags_setting}" = "xyes"],
		[
		LDFLAGS=""
		AX_APPEND_LINK_FLAGS([$ldflags_test], [LDFLAGS])
		LDFLAGS=$( echo $LDFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)



dnl ===============
dnl Check for C++11
dnl ===============
AX_CXX_COMPILE_STDCXX_11([], [mandatory])



dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.15 foreign dist-bzip2 no-dist-gzip subdir-objects silent-rules])



dnl =============
dnl Check for GSL
dnl =============
PKG_CHECK_MODULES([GSL], [gsl])



dnl ===============
dnl Check for Boost
dnl ===============
AX_BOOST_BASE([1.50], [], [AC_MSG_ERROR([could not find a suitable boost installation.])])
AX_BOOST_SYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_FILESYSTEM
AX_BOOST_THREAD



dnl ===============
dnl Check for Seqan
dnl ===============
AC_ARG_VAR([SEQAN_INCLUDEDIR], [contains the include path to the Seqan headers])
AC_MSG_CHECKING([whether seqan headers are available])
save_cppflags="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS -I$SEQAN_INCLUDEDIR"
AC_LINK_IFELSE([
	AC_LANG_SOURCE(
		[[
		#include <seqan/basic.h>
		int main() { seqan::Pair<int, int> pair; return 0; }
		]]
	)],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_ERROR([could not find seqan headers! If the seqan headers are not located within the compiler's default search path, specify them with the SEQAN_INCLUDEDIR variable to the configure script.])]
)
CPPFLAGS="$save_cppflags"



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	$PACKAGE_NAME $VERSION
	Using $ax_cv_cxx_compiler_vendor, $ax_cv_cxx_compiler_version

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
	CPPFLAGS:   $CPPFLAGS
	LDFLAGS:    $LDFLAGS
])

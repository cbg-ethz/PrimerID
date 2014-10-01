dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[PIDalyser],
	[0.1],
	[david.seifert@bsse.ethz.ch],
	[pidalyser],
	[http://www.cbg.ethz.ch/software/pidalyser])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG([C++])



dnl ========================================================================
dnl Check whether we want to set defaults for CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ========================================================================
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

AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [Disables setting of NDEBUG, which affects debugging asserts in Eigen])]
)




dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
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
		AX_APPEND_COMPILE_FLAGS([-O3 -Wall -pedantic], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_cppflags_setting}" = "xyes"],
		[
			CPPFLAGS=""
			AS_IF(
				[test "x$enable_debug" != "xyes"],
					[
					AX_APPEND_COMPILE_FLAGS([-DNDEBUG], [CPPFLAGS])
					]
			)
			CPPFLAGS=$( echo $CPPFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_ldflags_setting}" = "xyes"],
		[
		LDFLAGS=""
		AX_APPEND_LINK_FLAGS([-Wl,-O1 -Wl,--as-needed -Wl,-dead_strip_dylibs], [LDFLAGS])
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
AM_INIT_AUTOMAKE([1.14 foreign dist-bzip2 no-dist-gzip])
AM_MAINTAINER_MODE([enable])



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	$PACKAGE_NAME $VERSION

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
	CPPFLAGS:   $CPPFLAGS
	LDFLAGS:    $LDFLAGS
])

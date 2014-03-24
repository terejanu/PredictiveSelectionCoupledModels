# SYNOPSIS
#
#   Test for  
#
#   AM_PATH_QUESO([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-queso=DIR option. Searches --with-queso,
#   $QUESO_DIR, and the usual places for QUESO headers and libraries.
#
#   On success, sets QUESO_CFLAGS, QUESO_LIBS, and
#   #defines HAVE_QUESO.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2011-02-08 by Gabriel Terejanu
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_QUESO_ANN],
[

AC_ARG_VAR(QUESO_DIR,[root directory of QUESO installation])

AC_ARG_WITH(queso, 
  [AS_HELP_STRING([--with-queso[=DIR]],[root directory of QUESO installation (default = QUESO_DIR)])],
  [with_queso=$withval
if test "${with_queso}" != yes; then
    QUESO_PREFIX=$withval
fi
],[
with_queso=$withval
if test "x${QUESO_DIR}" != "x"; then
   QUESO_PREFIX=${QUESO_DIR}
fi
])

if test "${with_queso}" != no ; then

    if test -d "${QUESO_PREFIX}/lib" ; then
       QUESO_LIBS="-L${QUESO_PREFIX}/lib -lqueso -lANN"
    fi

    if test -d "${QUESO_PREFIX}/include" ; then
        dnl FIXME: CFLAGS should be reserved for C-compiler flags
        dnl FIXME: conflicts with use of CFLAGS below.
        QUESO_CFLAGS="-I${QUESO_PREFIX}/include"

	dnl PB: Added QUESO_CPPFLAGS
        QUESO_CPPFLAGS="-I${QUESO_PREFIX}/include"
    fi

    ac_QUESO_save_CFLAGS="$CFLAGS"
    ac_QUESO_save_CPPFLAGS="$CPPFLAGS"
    ac_QUESO_save_LDFLAGS="$LDFLAGS"
    ac_QUESO_save_LIBS="$LIBS"

    CFLAGS="${QUESO_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${QUESO_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${QUESO_LIBS} ${LDFLAGS}"


    AC_LANG_PUSH([C++])

    #--------------------------
    # Header check QUESO
    #--------------------------
    AC_CHECK_HEADER([queso.h], [found_header_QUESO=yes], [found_header_QUESO=no])

    #--------------------------
    # QUESO Library availability
    #--------------------------
    AC_MSG_CHECKING([for -lqueso linkage])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <uqEnvironment.h>],[QUESO::QUESO_get_numeric_version])],
    [TEST_LIBS="$TEST_LIBS -lqueso"] [
    AC_MSG_RESULT(yes)
    found_library_QUESO=yes ],[AC_MSG_RESULT(no)])

    #--------------------------
    # Header check ANN
    #--------------------------
    AC_CHECK_HEADER([ANN/ANN.h], [found_header_ANN=yes], [found_header_ANN=no])

    #--------------------------
    # Library availability
    #--------------------------
    AC_MSG_CHECKING([for -lANN linkage])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <ANN/ANN.h>],[annClose])],
    [TEST_LIBS="$TEST_LIBS -lANN"] [
    AC_MSG_RESULT(yes)
    found_library_ANN=yes ],[AC_MSG_RESULT(no)])

    AC_LANG_POP([C++])


    CFLAGS="$ac_QUESO_save_CFLAGS"
    CPPFLAGS="$ac_QUESO_save_CPPFLAGS"
    LDFLAGS="$ac_QUESO_save_LDFLAGS"
    LIBS="$ac_QUESO_save_LIBS"

    succeeded_QUESO=no
    if test "$found_header_QUESO" = yes; then
        if test "$found_library_QUESO" = yes; then
            succeeded_QUESO=yes
        fi
    fi

    succeeded_ANN=no
    if test "$found_header_ANN" = yes; then
        if test "$found_library_ANN" = yes; then
            succeeded_ANN=yes
        fi
    fi    

    if test "$succeeded_QUESO" = no; then
        ifelse([$2],,AC_MSG_ERROR([QUESO not found. Try either --with-queso or setting QUESO_DIR.]),
            [$2])
    else
        if test "$succeeded_ANN" = no; then
            ifelse([$2],,AC_MSG_ERROR([ANN not found. Make sure that QUESO is configured with --enable-ann=yes.]),
            [$2])
        else
            AC_DEFINE(HAVE_QUESO,1,[Define if QUESO is available])
            AC_SUBST(QUESO_CFLAGS)
            AC_SUBST(QUESO_CPPFLAGS)
            AC_SUBST(QUESO_LIBS)
            AC_SUBST(QUESO_PREFIX)
            ifelse([$1],,,[$1])
        fi
    fi

fi

])

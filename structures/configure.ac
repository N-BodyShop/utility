#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ(2.57)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([TipsyFile.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h endian.h machine/endian.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPE(u_int64_t, uint64_t)

# Checks for library functions.
AC_HEADER_STDC
AC_CHECK_FUNCS([floor memset sqrt xdr_hyper])

AC_OUTPUT

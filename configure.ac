dnl This file is part of epieov
dnl
dnl Copyright (C) 2020: University of Southern California
dnl                     Andrew D. Smith and Xiaojing Ji
dnl
dnl Authors: Andrew D. Smith
dnl
dnl This is free software: you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This software is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.

AC_INIT([epievo], [0.1.0], [andrewds@usc.edu],
        [epievo], [https://github.com/smithlabcode/epievo])
dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

dnl Less output to the terminal. Use "--disable-silent-rules" when
dnl running ./configure to reverse this.
AM_SILENT_RULES([yes])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_PROG_RANLIB

dnl check for required libraries
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas])
AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl])

AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT

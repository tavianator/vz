dnl Copyright (C) 2010 Tavian Barnes <tavianator@gmail.com>
dnl
dnl This file is part of The vZ Build Suite.
dnl
dnl The vZ Build Suite is free software; you can redistribute it and/or
dnl modify it under the terms of the GNU General Public License as
dnl published by the Free Software Foundation; either version 3 of the
dnl License, or (at your option) any later version.
dnl
dnl The vZ Build Suite is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program.  If not, see <http://www.gnu.org/licenses/>.

dnl Initialization
AC_PREREQ(2.68)
AC_INIT(vZ, 0.0.0, tavianator@gmail.com, vz,
        [http://www.gitorious.org/vz])
AM_INIT_AUTOMAKE([parallel-tests color-tests])
AM_SILENT_RULES([yes])

dnl Programs
AC_PROG_CXX
AC_PROG_LIBTOOL

dnl Generate Makefiles
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT

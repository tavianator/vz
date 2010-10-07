###########################################################################
## Copyright (C) 2009-2010 Tavian Barnes <tavianator@gmail.com>          ##
##                                                                       ##
## This file is part of The Dimension Build Suite.                       ##
##                                                                       ##
## The Dimension Build Suite is free software; you can redistribute it   ##
## and/or modify it under the terms of the GNU General Public License as ##
## published by the Free Software Foundation; either version 3 of the    ##
## License, or (at your option) any later version.                       ##
##                                                                       ##
## The Dimension Build Suite is distributed in the hope that it will be  ##
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty   ##
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  ##
## General Public License for more details.                              ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program.  If not, see <http://www.gnu.org/licenses/>. ##
###########################################################################

INCLUDES = -I$(top_srcdir)/src

check_PROGRAMS = Euler-test                                                    \
                 Midpoint-test                                                 \
                 Heun-test                                                     \
                 RK4-test                                                      \
                 HE12-test                                                     \
                 BS23-test                                                     \
                 RKF45-test                                                    \
                 CK45-test
TESTS          = $(check_PROGRAMS)

Euler_test_SOURCES    = Euler.cpp
Midpoint_test_SOURCES = Midpoint.cpp
Heun_test_SOURCES     = Heun.cpp
RK4_test_SOURCES      = RK4.cpp
HE12_test_SOURCES     = HE12.cpp
BS23_test_SOURCES     = BS23.cpp
RKF45_test_SOURCES    = RKF45.cpp
CK45_test_SOURCES     = CK45.cpp

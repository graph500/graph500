# -*- Makefile -*-
# Copyright 2010, University of Illinois at Urbana-Champaign
# See COPYING for license.
CFLAGS = -g -std=c99 -Wall
LDLIBS = -lm
# OSX does not support MAP_ANON in POSIX mode and the timers for MacOSX are
# not available for Leopard
# Use HAVE_MACH_ABSOLUTE_TIME to get the Mac OSX Timer
# Use HAVE_ALLOCA_H to include alloca.h
# Use HAVE_MPI_INT64_T to use the MPI 2.2 Datatype for int64_t items
# Use HAVE_UNISTD_H to include unistd.h (for getopt definitions)
CPPFLAGS = -DHAVE_MACH_ABSOLUTE_TIME -DHAVE_ALLOCA_H -DHAVE_MPI_INT64_T -DHAVE_UNISTD_H
BUILD_MPI = Yes
BUILD_OPENMP = No

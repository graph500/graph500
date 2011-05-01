/* Copyright (C) 2009-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef MPI_WORKAROUNDS_H
#define MPI_WORKAROUNDS_H

#include <mpi.h> /* To get version number */

/* Settings for user modification ----------------------------------- */

#if MPI_VERSION > 2 || (MPI_VERSION == 2 && MPI_SUBVERSION >= 2)
/* Use standard definitions in MPI 2.2. */
#else
/* Backup version for non-2.2-compliant MPI implementations. */
#error "Your MPI implementation is not compliant with the MPI 2.2 standard; please comment out this line in mpi_workarounds.h and ensure that the type definitions after it are correct."

#define FIND_MPI_INT_TYPE(t) (sizeof(t) == sizeof(signed char) ? MPI_SIGNED_CHAR : sizeof(t) == sizeof(short) ? MPI_SHORT : sizeof(t) == sizeof(int) ? MPI_INT : sizeof(t) == sizeof(long) ? MPI_LONG : sizeof(t) == sizeof(long long) ? MPI_LONG_LONG : MPI_DATATYPE_NULL)
#define FIND_MPI_UINT_TYPE(t) (sizeof(t) == sizeof(unsigned char) ? MPI_UNSIGNED_CHAR : sizeof(t) == sizeof(unsigned short) ? MPI_UNSIGNED_SHORT : sizeof(t) == sizeof(unsigned int) ? MPI_UNSIGNED : sizeof(t) == sizeof(unsigned long) ? MPI_UNSIGNED_LONG : sizeof(t) == sizeof(unsigned long long) ? MPI_UNSIGNED_LONG_LONG : MPI_DATATYPE_NULL)

#define MPI_INT64_T FIND_MPI_INT_TYPE(int64_t)
#define MPI_UINT64_T FIND_MPI_UINT_TYPE(uint64_t)
#define MPI_UINT32_T FIND_MPI_UINT_TYPE(uint32_t)
#define MPI_UINT16_T FIND_MPI_UINT_TYPE(uint16_t)
#define MPI_AINT FIND_MPI_UINT_TYPE(MPI_Aint)
#endif

/* If you have problems with one-sided operations (such as crashes in the
 * result validator), you can switch to an emulation of MPI-2 one-sided
 * operations: */
#undef EMULATE_ONE_SIDED
/* #define EMULATE_ONE_SIDED */

/* End of user settings ----------------------------------- */

#endif /* MPI_WORKAROUNDS_H */

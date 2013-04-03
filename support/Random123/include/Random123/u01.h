/*
Copyright 2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _random123_u01_dot_h_
#define _random123_u01_dot_h_

#include "features/compilerfeatures.h"

/** @defgroup u01_closed_open_W_M The u01 conversion functions

    These functions convert unsigned W-bit integers to real values
    (float or double) between 0.0 and 1.0 with mantissas of M bits.
    There are 12 functions, corresponding to the following choices:

     -  W = 32 or 64
     -  M = 24 or 53 (float or double)
     -  open0 or closed0 : whether the output is open or closed at 0.0
     -  open1 or closed1 : whether the output is open or closed at 1.0 

    The W=64 M=24 cases are not implemented.  To obtain an M=24 float
    from a uint64_t, use a cast (possibly with right-shift and bitwise
    and) to convert some of the bits of the uint64_t to a uint32_t and
    then use u01_x_y_32_24.  Note that the 64-bit random integers
    produced by the Random123 library are random in "all the bits", so
    with a little extra effort you can obtain two floats this way --
    one from the high bits and one from the low bits of the 64-bit
    value.

    If the output is open at one end, then the extreme
    value (0.0 or 1.0) will never be returned.  Conversely, if the output
    is closed at one end, then the extreme value is a possible
    return value.
    
    On x86 hardware, especially on 32bit machines, the use of
    internal 80bit x87-style floating point may result in
    'bonus' precision, which may cause closed intervals to not
    be really closed, i.e. the conversions below might not
    convert UINT{32,64}_MAX to 1.0.  This sort of issue is
    likely to occur when storing the output of a u01_*_32_24
    function in a double, though one can imagine getting extra
    precision artifacts when going from 64_53 as well.  Other
    artifacts may exist on some GPU hardware.  The tests in
    kat_u01_main.h try to expose such issues, but caveat emptor.

    @{
    @cond HIDDEN_FROM_DOXYGEN
 */

/* Hex floats were standardized by C in 1999, but weren't standardized
   by C++ until 2011.  So, we're obliged to write out our constants in
   decimal, even though they're most naturally expressed in binary.
   We cross our fingers and hope that the compiler does the compile-time
   constant arithmetic properly.
*/
#define R123_0x1p_32f (1.f/4294967296.f)
#define R123_0x1p_24f (1.f/16777216.f)
#define R123_0x1fffffep_25f (16777215.f * R123_0x1p_24f * R123_0x1p_24f)
#define R123_0x1p_64 (1./(4294967296.*4294967296.))
#define R123_0x1p_53 (1./(4294967296.*2097152.))
#define R123_0x1fffffffffffffp_54 (9007199254740991.*R123_0x1p_53*R123_0x1p_53)
#define R123_0x1p_32 (1./4294967296.)
#define R123_0x100000001p_32 (4294967297.*R123_0x1p_32*R123_0x1p_32)

/** @endcond */

#ifdef __cplusplus
extern "C"{
#endif

/* narrowing conversions:  uint32_t to float */
R123_CUDA_DEVICE R123_STATIC_INLINE float u01_closed_closed_32_24(uint32_t i){
    return i*R123_0x1p_32f; /* 0x1.p-32f */
}

R123_CUDA_DEVICE R123_STATIC_INLINE float u01_closed_open_32_24(uint32_t i){
    return (i>>8)*R123_0x1p_24f; /* 0x1.0p-24f; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE float u01_open_closed_32_24(uint32_t i){
    return (1+(i>>8))*R123_0x1p_24f; /* *0x1.0p-24f; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE float u01_open_open_32_24(uint32_t i){
    return (1+(i>>8))*R123_0x1fffffep_25f; /* 0x1.fffffep-25f; */
}

#if R123_USE_U01_DOUBLE
/* narrowing conversions:  uint64_t to double */
R123_CUDA_DEVICE R123_STATIC_INLINE double u01_closed_closed_64_53(uint64_t i){
    return i*R123_0x1p_64; /* 0x1.p-64; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE double u01_closed_open_64_53(uint64_t i){
    return (i>>11)*R123_0x1p_53; /* 0x1.0p-53; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE double u01_open_closed_64_53(uint64_t i){
    return (1+(i>>11))*R123_0x1p_53; /* 0x1.0p-53; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE double u01_open_open_64_53(uint64_t i){
    return (1+(i>>11))*R123_0x1fffffffffffffp_54; /* 0x1.fffffffffffffp-54; */
}

/* widening conversions:  u32 to double */
R123_CUDA_DEVICE R123_STATIC_INLINE double u01_closed_closed_32_53(uint32_t i){
    return i*R123_0x100000001p_32; /* 0x1.00000001p-32; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE double u01_closed_open_32_53(uint32_t i){
    return i*R123_0x1p_32; /* 0x1.p-32; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE double u01_open_closed_32_53(uint32_t i){
    return (1.+i)*R123_0x1p_32; /* 0x1.p-32; */
}

R123_CUDA_DEVICE R123_STATIC_INLINE double u01_open_open_32_53(uint32_t i){
    return (0.5+i)*R123_0x1p_32; /* 0x1.p-32; */
}
#endif /* R123_USE_U01_DOUBLE */

#ifdef __cplusplus
}
#endif

/** @} */
#endif

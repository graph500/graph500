/*
Copyright 2012, Georgia Institute of Technology.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of Georgia Institute of Technology nor the names of
  its contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

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
#ifndef __xlcfeatures_dot_hpp
#define __xlcfeatures_dot_hpp

#ifndef R123_STATIC_INLINE
#define R123_STATIC_INLINE static __inline
#endif

/* Needs xlc version check, etc. */
#ifndef R123_FORCE_INLINE
#define R123_FORCE_INLINE(decl) decl __attribute__((always_inline))
#endif

#ifndef R123_CUDA_DEVICE
#define R123_CUDA_DEVICE
#endif

#ifndef R123_ASSERT
#include <assert.h>
#define R123_ASSERT(x) assert(x)
#endif

#ifndef R123_BUILTIN_EXPECT
#define R123_BUILTIN_EXPECT(expr,likely) __builtin_expect(expr,likely)
#endif

#ifndef R123_USE_CXX0X
#define R123_USE_CXX0X 0
#endif

#ifndef R123_USE_AES_NI
#define R123_USE_AES_NI 0
#endif

#ifndef R123_USE_SSE4_2
#define R123_USE_SSE4_2 0
#endif

#ifndef R123_USE_SSE4_1
#define R123_USE_SSE4_1 0
#endif

#ifndef R123_USE_SSE
#define R123_USE_SSE 0
#endif

#ifndef R123_USE_STD_RANDOM
#define R123_USE_STD_RANDOM 1
#endif

#ifndef R123_USE_AES_OPENSSL
#define R123_USE_AES_OPENSSL 0
#endif

#ifndef R123_USE_GNU_UINT128
#define R123_USE_GNU_UINT128 0
#endif

#ifndef R123_USE_ASM_GNU
#define R123_USE_ASM_GNU 0
#endif

#ifndef R123_USE_CPUID_MSVC
#define R123_USE_CPUID_MSVC 0
#endif

#ifndef R123_USE_X86INTRIN_H
#define R123_USE_X86INTRIN_H 0
#endif

#ifndef R123_USE_IA32INTRIN_H
#define R123_USE_IA32INTRIN_H 0
#endif

#ifndef R123_USE_EMMINTRIN_H
#define R123_USE_EMMINTRIN_H 0
#endif

#ifndef R123_USE_SMMINTRIN_H
#define R123_USE_SMMINTRIN_H 0
#endif

#ifndef R123_USE_WMMINTRIN_H
#define R123_USE_WMMINTRIN_H 0
#endif

#ifndef R123_USE_ALTIVEC_H
#define R123_USE_ALTIVEC_H 0
#endif

#ifndef R123_USE_INTRIN_H
#define R123_USE_INTRIN_H 0
#endif

#ifndef R123_USE_MULHILO16_ASM
#define R123_USE_MULHILO16_ASM 0
#endif

#ifndef R123_USE_MULHILO32_ASM
#define R123_USE_MULHILO32_ASM 0
#endif

#ifndef R123_USE_MULHILO64_ASM
#define R123_USE_MULHILO64_ASM 0
#endif

#ifndef R123_USE_MULHILO64_MSVC_INTRIN
#define R123_USE_MULHILO64_MSVC_INTRIN 0
#endif

#ifndef R123_USE_MULHILO64_CUDA_INTRIN
#define R123_USE_MULHILO64_CUDA_INTRIN 0
#endif

#ifndef R123_USE_MULHILO64_OPENCL_INTRIN
#define R123_USE_MULHILO64_OPENCL_INTRIN 0
#endif

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include <stdint.h>
#ifndef UINT64_C
#error UINT64_C not defined.  You must define __STDC_CONSTANT_MACROS before you #include <stdint.h>
#endif

// If you add something, it must go in all the other XXfeatures.hpp
// and in ../ut_features.cpp
#endif

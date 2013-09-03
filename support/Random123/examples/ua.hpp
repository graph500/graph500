/*
Copyright 2010-2011, D. E. Shaw Research.
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

#if __cplusplus<201103L
#error "This file uses C++11 features: <array> and <type_traits> and trailing return types"
#endif

#ifndef __r123_ua_dot_hpp
#define __r123_ua_dot_hpp
// uniform.hpp provides functions for converting single integers
// to uniformly distributed reals.
//
// This file takes r123arrays of integers and converts them into
// std::arrays of floating point type.  Why std::array and not
// r123array?  Because with r123array, the size is baked into the name
// and is not a template parameter.  Since it's not a template
// parameter, we can't manufacture an r123array of the right size
// during template instantiation.
//
// This file may not be as portable, and has not been tested as
// rigorously as the files in the library itself, i.e., those in
// ../include/Random123.  Nevertheless, we hope it is useful and we
// encourage developers to copy it and modify it for their own
// use.  We invite comments and improvements.

#include <Random123/features/compilerfeatures.h>
#include <array>
#include "uniform.hpp"

namespace r123{

template <typename Ftype, typename IAtype>
static inline auto ua01(IAtype in) -> 
    std::array<Ftype, IAtype::static_size>
{
    std::array<Ftype, IAtype::static_size> ret;
    size_t i=0;
    for(auto e : in){
        ret[i++] = u01<Ftype>(e);
    }
    return ret;
}

template <typename Ftype, typename IAtype>
static inline auto uaneg11(IAtype in) -> 
    std::array<Ftype, IAtype::static_size>
{
    std::array<Ftype, IAtype::static_size> ret;
    size_t i=0;
    for(auto e : in){
        ret[i++] = uneg11<Ftype>(e);
    }
    return ret;
}

template <typename Ftype, typename IAtype>
static inline auto ua01fixedpt(IAtype in) -> 
    std::array<Ftype, IAtype::static_size>
{
    std::array<Ftype, IAtype::static_size> ret;
    size_t i=0;
    for(auto e : in){
        ret[i++] = u01fixedpt<Ftype>(e);
    }
    return ret;
}


} // namespace r123

#endif


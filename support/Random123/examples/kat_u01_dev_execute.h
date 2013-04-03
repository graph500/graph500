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
#include "kat_u01.h"

#ifndef KAT_KERNEL
#define KAT_KERNEL
#endif

#ifndef KAT_GLOBAL
#define KAT_GLOBAL
#endif

#ifndef KAT_THREADID
#define KAT_THREADID    (0)
#endif

/*
 * convert to floats and doubles
 */
KAT_KERNEL void
dev_execute_tests(KAT_GLOBAL uint64_t *vals, KAT_GLOBAL KatU01Result *ret) {
    unsigned tid = KAT_THREADID;
    uint64_t v64 = vals[tid];
    uint32_t v32 = vals[tid] & 0xffffffff;

    if (v32 == v64) {
	ret[tid].valuef[0] = u01_closed_closed_32_24(v32);
	ret[tid].valuef[1] = u01_closed_open_32_24(v32);
	ret[tid].valuef[2] = u01_open_closed_32_24(v32);
	ret[tid].valuef[3] = u01_open_open_32_24(v32);

	ret[tid].valued[0] = u01_closed_closed_32_53(v32);
	ret[tid].valued[1] = u01_closed_open_32_53(v32);
	ret[tid].valued[2] = u01_open_closed_32_53(v32);
	ret[tid].valued[3] = u01_open_open_32_53(v32);
    } else {
	int i;
	for (i = 0; i < 4; i++) {
	    ret[tid].valuef[i] = 0.f;
	    ret[tid].valued[i] = 0.;
	}
    }

    ret[tid].valued[4] = u01_closed_closed_64_53(v64);
    ret[tid].valued[5] = u01_closed_open_64_53(v64);
    ret[tid].valued[6] = u01_open_closed_64_53(v64);
    ret[tid].valued[7] = u01_open_open_64_53(v64);
}

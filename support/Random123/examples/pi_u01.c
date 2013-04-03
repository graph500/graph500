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
#include <stdio.h>
#include "Random123/threefry.h"
#include "Random123/u01.h"

/* Compute pi, using the u01 conversion with threefry2x64 and threefry2x32 */

#include "pi_check.h"

int main(int argc, char **argv){
    unsigned long hits = 0, tries = 0;
    int errs = 0;

    threefry2x64_ctr_t c = {{0}}, r;
    threefry2x64_key_t k = {{R123_64BIT(0xdeadbeef12345678)}};
    threefry2x32_ctr_t ch = {{0}}, rh;
    threefry2x32_key_t kh = {{0xdecafbad}};
    (void)argc; (void)argv; /* unused */
    printf("%lu uniform doubles from threefry2x64\n", NTRIES);
    while (tries < NTRIES) {
            double x, y;
            c.v[0]++; /* increment the counter */
	    r = threefry2x64(c, k);
            x = 2.*u01_open_open_64_53(r.v[0]) - 1.;
            y = 2.*u01_open_open_64_53(r.v[1]) - 1.;
            if( x*x + y*y < 1.0 )
                hits++;
	    tries++;
    }
    errs += pi_check(hits, tries);

    printf("%lu uniform doubles from threefry2x32\n", NTRIES);
    hits = tries = 0;
    while (tries < NTRIES) {
            double x, y;
            ch.v[0]++; /* increment the counter */
	    rh = threefry2x32(ch, kh);
            x = 2.*u01_open_open_32_53(rh.v[0]) - 1.;
            y = 2.*u01_open_open_32_53(rh.v[1]) - 1.;
            if( x*x + y*y < 1.0 )
                hits++;
	    tries++;
    }
    errs += pi_check(hits, tries);

    printf("%lu uniform floats from threefry2x32\n", NTRIES);
    hits = tries = 0;
    while (tries < NTRIES) {
            float x, y;
            ch.v[0]++; /* increment the counter */
	    rh = threefry2x32(ch, kh);
            x = 2.f*u01_open_open_32_24(rh.v[0]) - 1.f;
            y = 2.f*u01_open_open_32_24(rh.v[1]) - 1.f;
            if( x*x + y*y < 1.0 )
                hits++;
	    tries++;
    }
    errs += pi_check(hits, tries);
    return errs;
}

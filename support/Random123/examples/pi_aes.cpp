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
#include <Random123/aes.h>
#include <stdio.h>
#include <assert.h>
using namespace r123;

// Everyone's favorite PRNG example: calculate pi/4 by throwing darts
// at a square board and counting the fraction that are inside the
// inscribed circle.

// This version uses the C++ API to AESNI.

#include "pi_check.h"

int main(int, char **){
#if R123_USE_AES_NI
    unsigned long hits = 0, tries = 0;
    const int64_t two_to_the_62 = ((int64_t)1)<<62;

    if (!haveAESNI()) {
	std::cerr << "AES-NI instructions not available on this hardware, skipping the pi_aes test." << std::endl;
        return 0;
    }
    typedef AESNI4x32 G;
    G generator;
    G::ukey_type ukey = {{0x11111111, 0x22222222, 0x33333333, 0x44444444}};
    // The key_type constructor transforms the 128bit AES ukey_type to an expanded (1408bit) form.
    G::key_type key = ukey;
    G::ctr_type ctr = {{0xdeadbeef, 0xbeadcafe, 0x12345678, 0x90abcdef}};

    printf("Throwing %lu darts at a square board using AESNI4x32\n", NTRIES);
    std::cout << "Initializing AES key with userkey: " << std::hex << ukey << " ctr: " << ctr << std::endl;

    while(tries < NTRIES){
        ctr.incr();
        G::ctr_type r = generator(ctr, key);
	if (tries == 0) {
	    std::cout << "first random from AESNI is " << std::hex << r << std::endl;;
	}
        for(size_t j=0; j<r.size(); j+=2){
            int64_t x = (int32_t)r[j];
            int64_t y = (int32_t)r[j+1];
            if( (x*x + y*y) < two_to_the_62 )
                hits++;
            tries++;
        }
    }
    return pi_check(hits, tries);
#else
    std::cout << "AESNI RNG not compiled into this binary, skipping the pi_aes test.\n";
    return 0;	// Not a failure to not have AESNI compiled into this.
#endif
}

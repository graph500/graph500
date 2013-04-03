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

#if 0 /* change this to 1 to test the failure code paths by adding a bad test */
#define KAT_U01_TEST_FAIL
#endif

/*
 * Define to %a if your compiler vendor has gotten around to
 * supporting %a.  Congratulations!
 */
#ifndef KAT_U01_FMT
#define KAT_U01_FMT "%.17f"
#endif

#include "util.h"
#include "kat_u01.h"

/*
 * The order of expected is {c-c,c-o,o-c,o-o}{32_{24,53},64_53},
 * see dev_execute_tests() in kat_u01_dev_execute.c.
 * The tests below were generated with
 * run_u01_c.sh 0 1 0xff 0x100 0x319a58f11e1dc63c 0xc611918ed445030f 0x12d8424b 0xa3b99ed1
 */
KatU01Test tests[] = {
#ifdef KAT_U01_TEST_FAIL
/* deliberately incorrect test outputs */
{R123_64BIT(1), {"0", "1", "0", "0",
	"0", "0", "0", "0",
	"0", "1", "0", "0"}},
#endif
{ R123_64BIT(0), { "0x0p+0","0x0p+0","0x1p-24","0x1.fffffep-25",
	"0x0p+0","0x0p+0","0x1p-32","0x1p-33",
	"0x0p+0","0x0p+0","0x1p-53","0x1.fffffffffffffp-54",}},
{ R123_64BIT(0xffffffff), { "0x1p+0","0x1.fffffep-1","0x1p+0","0x1.fffffep-1",
	"0x1p+0","0x1.fffffffep-1","0x1p+0","0x1.ffffffffp-1",
	"0x1.fffffffep-33","0x1.fffffp-33","0x1p-32","0x1.fffffffffffffp-33",}},
{ R123_64BIT(0xffffffffffffffff), { "","","","",
	"","","","",
	"0x1p+0","0x1.fffffffffffffp-1","0x1p+0","0x1.fffffffffffffp-1",}},
{ R123_64BIT(0x1), { "0x1p-32","0x0p+0","0x1p-24","0x1.fffffep-25",
	"0x1.00000001p-32","0x1p-32","0x1p-31","0x1.8p-32",
	"0x1p-64","0x0p+0","0x1p-53","0x1.fffffffffffffp-54",}},
{ R123_64BIT(0xfffffffe), { "0x1p+0","0x1.fffffep-1","0x1p+0","0x1.fffffep-1",
	"0x1.fffffffep-1","0x1.fffffffcp-1","0x1.fffffffep-1","0x1.fffffffdp-1",
	"0x1.fffffffcp-33","0x1.fffffp-33","0x1p-32","0x1.fffffffffffffp-33",}},
{ R123_64BIT(0xfffffffffffffffe), { "","","","",
	"","","","",
	"0x1p+0","0x1.fffffffffffffp-1","0x1p+0","0x1.fffffffffffffp-1",}},
{ R123_64BIT(0xff), { "0x1.fep-25","0x0p+0","0x1p-24","0x1.fffffep-25",
	"0x1.fe000001fep-25","0x1.fep-25","0x1p-24","0x1.ffp-25",
	"0x1.fep-57","0x0p+0","0x1p-53","0x1.fffffffffffffp-54",}},
{ R123_64BIT(0xffffff00), { "0x1.fffffep-1","0x1.fffffep-1","0x1p+0","0x1.fffffep-1",
	"0x1.fffffe02p-1","0x1.fffffep-1","0x1.fffffe02p-1","0x1.fffffe01p-1",
	"0x1.fffffep-33","0x1.fffffp-33","0x1p-32","0x1.fffffffffffffp-33",}},
{ R123_64BIT(0xffffffffffffff00), { "","","","",
	"","","","",
	"0x1p+0","0x1.fffffffffffffp-1","0x1p+0","0x1.fffffffffffffp-1",}},
{ R123_64BIT(0x100), { "0x1p-24","0x1p-24","0x1p-23","0x1.fffffep-24",
	"0x1.00000001p-24","0x1p-24","0x1.01p-24","0x1.008p-24",
	"0x1p-56","0x0p+0","0x1p-53","0x1.fffffffffffffp-54",}},
{ R123_64BIT(0xfffffeff), { "0x1.fffffep-1","0x1.fffffcp-1","0x1.fffffep-1","0x1.fffffcp-1",
	"0x1.fffffep-1","0x1.fffffdfep-1","0x1.fffffep-1","0x1.fffffdffp-1",
	"0x1.fffffdfep-33","0x1.fffffp-33","0x1p-32","0x1.fffffffffffffp-33",}},
{ R123_64BIT(0xfffffffffffffeff), { "","","","",
	"","","","",
	"0x1p+0","0x1.fffffffffffffp-1","0x1p+0","0x1.fffffffffffffp-1",}},
{ R123_64BIT(0x319a58f11e1dc63c), { "","","","",
	"","","","",
	"0x1.8cd2c788f0ee3p-3","0x1.8cd2c788f0eep-3","0x1.8cd2c788f0ee4p-3","0x1.8cd2c788f0ee3p-3",}},
{ R123_64BIT(0xce65a70ee1e239c3), { "","","","",
	"","","","",
	"0x1.9ccb4e1dc3c47p-1","0x1.9ccb4e1dc3c47p-1","0x1.9ccb4e1dc3c48p-1","0x1.9ccb4e1dc3c47p-1",}},
{ R123_64BIT(0xc611918ed445030f), { "","","","",
	"","","","",
	"0x1.8c23231da88ap-1","0x1.8c23231da88ap-1","0x1.8c23231da88a1p-1","0x1.8c23231da88ap-1",}},
{ R123_64BIT(0x39ee6e712bbafcf0), { "","","","",
	"","","","",
	"0x1.cf7373895dd7ep-3","0x1.cf7373895dd7cp-3","0x1.cf7373895dd8p-3","0x1.cf7373895dd7fp-3",}},
{ R123_64BIT(0x12d8424b), { "0x1.2d8424p-4","0x1.2d842p-4","0x1.2d843p-4","0x1.2d842ep-4",
	"0x1.2d8424b12d842p-4","0x1.2d8424bp-4","0x1.2d8424cp-4","0x1.2d8424b8p-4",
	"0x1.2d8424bp-36","0x1.2d84p-36","0x1.2d848p-36","0x1.2d847ffffffffp-36",}},
{ R123_64BIT(0xed27bdb4), { "0x1.da4f7cp-1","0x1.da4f7ap-1","0x1.da4f7cp-1","0x1.da4f7ap-1",
	"0x1.da4f7b69da4f8p-1","0x1.da4f7b68p-1","0x1.da4f7b6ap-1","0x1.da4f7b69p-1",
	"0x1.da4f7b68p-33","0x1.da4f7p-33","0x1.da4f8p-33","0x1.da4f7ffffffffp-33",}},
{ R123_64BIT(0xffffffffed27bdb4), { "","","","",
	"","","","",
	"0x1.ffffffffda4f8p-1","0x1.ffffffffda4f7p-1","0x1.ffffffffda4f8p-1","0x1.ffffffffda4f7p-1",}},
{ R123_64BIT(0xa3b99ed1), { "0x1.47733ep-1","0x1.47733cp-1","0x1.47733ep-1","0x1.47733cp-1",
	"0x1.47733da347734p-1","0x1.47733da2p-1","0x1.47733da4p-1","0x1.47733da3p-1",
	"0x1.47733da2p-33","0x1.47733p-33","0x1.47734p-33","0x1.47733ffffffffp-33",}},
{ R123_64BIT(0x5c46612e), { "0x1.711984p-2","0x1.711984p-2","0x1.711988p-2","0x1.711986p-2",
	"0x1.711984b971198p-2","0x1.711984b8p-2","0x1.711984bcp-2","0x1.711984bap-2",
	"0x1.711984b8p-34","0x1.71198p-34","0x1.7119ap-34","0x1.71199ffffffffp-34",}},
{ R123_64BIT(0xffffffff5c46612e), { "","","","",
	"","","","",
	"0x1.fffffffeb88ccp-1","0x1.fffffffeb88ccp-1","0x1.fffffffeb88cdp-1","0x1.fffffffeb88ccp-1",}},
};

int check_results(KatU01Test *tests, size_t ntests, KatU01Result *results)
{
    double d;
    float f;
    size_t i, j;
    int errs = 0;
    uint32_t inp32;
    uint64_t inp;

    for (i = 0; i < ntests; i++) {
	if (debug) {
	    printf("tests[%lu] = {\n", (unsigned long)i);
	    for (j = 0; j < NUM_CONVERTERS; j++) {
		printf("\"%s\" " KAT_U01_FMT " , ", tests[i].expected[j], hextod(tests[i].expected[j]));
	    }
	    printf("}\n");
	}
	inp = tests[i].input;
	inp32 = inp & 0xffffffff;
	if (inp32 == inp) {
	    /*
	     * On IA-32, it is important to store the results from
	     * u01*32_24 in float if you want exact matches to
	     * what one would expect from IEEE math.  If one
	     * stores them in double, then it's likely that one
	     * will get more bits than expected, thanks to 80bit
	     * internal numbers.
	     */
	    for (j = 0; j < NUM_FLOAT_CONVERTERS; j++) {
		f = hextod(tests[i].expected[j]);
		dprintf(("valuef %lu %lu " KAT_U01_FMT " " KAT_U01_FMT "\n", (unsigned long)i, (unsigned long)j, results[i].valuef[j], f));
		if (f != results[i].valuef[j]) {
		    fprintf(stderr, "*** error test %lu valuef[%lu] got " KAT_U01_FMT " expected " KAT_U01_FMT "\n", (unsigned long)i, (unsigned long)j, results[i].valuef[j], f);
		    errs++;
		}
	    }
	}
	for (j = (inp32 == inp) ? 0 : 4; j < NUM_DOUBLE_CONVERTERS; j++) {
	    d = hextod(tests[i].expected[j + NUM_FLOAT_CONVERTERS]);
	    dprintf(("valued %lu %lu " KAT_U01_FMT " " KAT_U01_FMT "\n", (unsigned long)i, (unsigned long)j, results[i].valued[j], d));
	    if (d != results[i].valued[j]) {
		fprintf(stderr, "*** error test %lu valued[%lu] got " KAT_U01_FMT " expected " KAT_U01_FMT "\n", (unsigned long)i, (unsigned long)j, results[i].valued[j], d);
		errs++;
	    }
	}
    }
    return errs;
}

extern void host_execute_tests(uint64_t *, size_t, KatU01Result *);

int debug = 0;
int verbose = 0;
const char *progname;

int
main(int argc, char **argv)
{
    char *cp;
    int errs = 0;
    size_t i, n = 0;
    uint64_t *inp;
    KatU01Result *out;
    
    (void)argc; /* unused */
    progname = argv[0];
    if ((cp = getenv("KAT_U01_DEBUG")) != NULL) debug = atoi(cp);
    /* First test that hextod appears to work */
    assert(hextod("0") == 0.);
    assert(hextod("0x0p0") == 0.);
    assert(hextod("0x0.0p0") == 0.);
    assert(hextod("1") == 1.);
    assert(hextod("0x1p+0") == 1.);
    assert(hextod("0x1.0p+0") == 1.);
    assert(hextod("0xe.1p-4") == 0.87890625);
    assert(hextod("0xf.ap+6") == 1000.);
    assert(hextod("0xf.fffffffp+28") == 4294967295.);
    assert(hextod("0xb.fffffffffffffcp+54") == 216172782113783807.);

    n = sizeof(tests)/sizeof(tests[0]);
#ifdef KAT_U01_TEST_FAIL
    printf("KAT_U01_TEST_FAIL defined, expect errors for test 0, valuef[0..%d] and valued[0..%d], i.e. %d test values will fail\n",
	   NUM_FLOAT_CONVERTERS-1, NUM_DOUBLE_CONVERTERS-1, NUM_CONVERTERS);
    fflush(stdout);
#endif
    CHECKNOTZERO(inp = (uint64_t *)malloc(n*sizeof(*inp)));
    for (i = 0; i < n; i++) {
	inp[i] = tests[i].input;
    }
    CHECKNOTZERO(out = (KatU01Result *)malloc(n*sizeof(*out)));
    host_execute_tests(inp, n, out);
    errs = check_results(tests, n, out);
    if (errs == 0) {
	printf("%lu test values PASSED\n", (unsigned long)n*NUM_CONVERTERS);
    } else {
	printf("%d test values FAILED (out of %lu)\n", errs, (unsigned long)n*NUM_CONVERTERS);
	return 1;
    }
    free(inp);
    free(out);
    return 0;
}


#include <inttypes.h>
#include <assert.h>

#define IN_GLOBALS_C
#include "globals.h"

int SCALE = 0, EF, NROOT, MAXWEIGHT;
int64_t NV, NE, Z, Zinv;
float A, B, NOISEFACT;

static int64_t
extended_gcd (int64_t a, int64_t b,
	      int64_t * restrict x_out, int64_t * restrict y_out)
{
  assert (a > 0);
  assert (b > 0);

  int64_t x = 0, y = 1, prevx = 1, prevy = 0;

  while (b) {
    const int64_t q = a/b;
    const int64_t r = a%b;
    int64_t t;

    a = b;
    b = r;
    t = prevx - q * x;
    prevx = x;
    x = t;
    t = prevy - q * y;
    prevy = y;
    y = t;
  }
  *x_out = prevx;
  *y_out = prevy;
  return a;
}

void
init_globals (int scale, int ef, int maxweight, int nroot,
	      float a, float b, float noisefact)
{
  SCALE = scale;
  NV = ((int64_t)1) << SCALE;
  EF = ef;
  NE = NV * ef;
  MAXWEIGHT = maxweight;
  NROOT = (nroot > NV? NV : (nroot > NROOT_MAX? NROOT_MAX : nroot));
  Z = Zinv = -1;
  A = a;
  B = b;
  NOISEFACT = noisefact;
  for (int64_t k = (3*NE)/4; k < NE; ++k) {
    int64_t gcd, x, y;
    gcd = extended_gcd (k, NE, &x, &y);
    if (1 == gcd) {
      Z = k;
      Zinv = x;
      break;
    }
  }
  assert (Z > 0);
  assert (Zinv > 0);
  assert (1 == (Z*Zinv) % NE);
}


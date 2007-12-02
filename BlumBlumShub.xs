#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <gmp.h>

void bbs(mpz_t * outref, mpz_t * p, mpz_t * q, mpz_t * seed, int bits_required) {
     mpz_t n, gcd, one;
     unsigned long i, k;

     if(mpz_fdiv_ui(*p, 4) != 3) croak ("First prime is unsuitable for Blum-Blum-Shub prbg (must be congruent to 3, mod 4)");
     if(mpz_fdiv_ui(*q, 4) != 3) croak ("Second prime is unsuitable for Blum-Blum-Shub prbg (must be congruent to 3, mod 4)"); 
     mpz_init(n);

     mpz_mul(n, *p, *q);

     if(mpz_cmp_ui(*seed, 0) < 0)croak("Negative seed supplied to bbs function");
     if(mpz_cmp(*seed, n) >= 0)croak("Seed supplied to bbs function is too big");

     mpz_init(gcd);

     mpz_gcd(gcd, *seed, n);
     if(mpz_cmp_ui(gcd, 1)) croak("gcd(seed, p * q) != 1");

     mpz_powm_ui(*seed, *seed, 2, n);

     mpz_init_set_ui(*outref, 0);
     mpz_init_set_ui(one, 1);

     for(i = 0; i < bits_required; ++i) {
         mpz_powm_ui(*seed, *seed, 2, n);
         k = mpz_tstbit(*seed, 0);
         if(k) {
            mpz_mul_2exp(gcd, one, i);      
            mpz_add(*outref, gcd, *outref);
            }
         }

     mpz_clear(n); 
     mpz_clear(gcd);
     mpz_clear(one);

}

void bbs_seedgen(mpz_t * seed, mpz_t * p, mpz_t * q) {
     mpz_t n, gcd;
     gmp_randstate_t state;

     mpz_init(n);

     mpz_mul(n, *p, *q);

     mpz_init(gcd);

     if(mpz_cmp_ui(*seed, 0) < 0)croak("Negative seed supplied to bbs_seedgen");
     gmp_randinit_default(state);
     gmp_randseed(state, *seed);
     mpz_urandomm(*seed, state, n);
     gmp_randclear(state);

     while(1) {
       mpz_gcd(gcd, *seed, n);
       if(!mpz_cmp_ui(gcd, 1)) break;
       mpz_sub_ui(*seed, *seed, 1);
     }

     mpz_clear(n);
     mpz_clear(gcd);
}

int monobit(mpz_t * bitstream) {
    unsigned long len, i, count = 0;

    len = mpz_sizeinbase(*bitstream, 2);

    if(len > 20000) croak("Wrong size random sequence for monobit test");
    if(len < 19967) {
       warn("More than 33 leading zeroes in monobit test\n");
       return 0;
       }

    count = mpz_popcount(*bitstream);

    if(count > 9654 && count < 10346) return 1;
    return 0;

}

int longrun(mpz_t * bitstream) {
    unsigned long i, el, init = 0, count = 0, len, t;

    len = mpz_sizeinbase(*bitstream, 2);

    if(len > 20000) croak("Wrong size random sequence for longrun test");
    if(len < 19967) {
       warn("More than 33 leading zeroes in longrun test\n");
       return 0;
       }

    el = mpz_tstbit(*bitstream, 0);

    for(i = 0; i < len; ++i) {
        t = mpz_tstbit(*bitstream, i);
        if(t == el) ++count;
        else {
           el = t;
           if(count > init) init = count;
           count = 1;
           }
        }

    if(init < 34 && count < 34) return 1;
    return 0;

}

int runs(mpz_t * bitstream) {
    int b[6] = {0,0,0,0,0,0}, g[6] = {0,0,0,0,0,0},
    len, count = 1, i, t, diff;

    len = mpz_sizeinbase(*bitstream, 2);
    diff = 20000 - len;

    if(len > 20000) croak("Wrong size random sequence for runs test");
    if(len < 19967) {
       warn("More than 33 leading zeroes in runs test\n");
       return 0;
       }

    --len;

    for(i = 0; i < len; ++i) {
        t = mpz_tstbit(*bitstream, i);
        if(t == mpz_tstbit(*bitstream, i + 1)) ++ count;
        else {
           if(t) {
              if(count >= 6) ++b[5];
              else ++b[count - 1];
              }
            else {
              if(count >= 6) ++g[5];
              else ++g[count - 1];
              }
            count = 1;
            }
         }

     if(count >= 6) {
        if(mpz_tstbit(*bitstream, len)) {
           ++b[5];
           if(diff) ++g[diff - 1];
           }
        else ++g[5];
        }
     else {
        if(mpz_tstbit(*bitstream, len)) {
           ++b[count - 1];
           if(diff) ++g[diff - 1];
           }
        else {
          count += diff;
          if(count >= 6) ++g[5];
          else ++g[count - 1];
          }
        }

             
        if (
            b[0] <= 2267 || g[0] <= 2267 ||
            b[0] >= 2733 || g[0] >= 2733 ||
            b[1] <= 1079 || g[1] <= 1079 ||
            b[1] >= 1421 || g[1] >= 1421 ||
            b[2] <= 502  || g[2] <= 502  ||
            b[2] >= 748  || g[2] >= 748  ||
            b[3] <= 223  || g[3] <= 223  ||
            b[3] >= 402  || g[3] >= 402  ||
            b[4] <= 90   || g[4] <= 90   ||
            b[4] >= 223  || g[4] >= 223  ||
            b[5] <= 90   || g[5] <= 90   ||
            b[5] >= 223  || g[5] >= 223
           ) return 0;

    return 1;

}

int poker (mpz_t * bitstream) {
    int counts[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
       len, i, st, diff;
    double n = 0;

    len = mpz_sizeinbase(*bitstream, 2);

    if(len > 20000) croak("Wrong size random sequence for poker test");
    if(len < 19967) {
       warn("More than 33 leading zeroes in poker test\n");
       return 0;
       }

/* pad with trailing zeroes (if necessary) to achieve length of 20000 bits. */
    diff = 20000 - len;
    if(diff) mpz_mul_2exp(*bitstream, *bitstream, diff);
    if(mpz_sizeinbase(*bitstream, 2) != 20000) croak("Bug in bit sequence manipulation in poker() function");

    for(i = 0; i < 20000; i += 4) {
        st = mpz_tstbit(*bitstream, i) +
             (mpz_tstbit(*bitstream, i + 1) * 2) +
             (mpz_tstbit(*bitstream, i + 2) * 4) +
             (mpz_tstbit(*bitstream, i + 3) * 8);
        ++counts[st];
        } 


    for(i = 0; i < 16; ++i) n += counts[i] * counts[i];

    n /= 5000;
    n *= 16;
    n -= 5000;

    if(n > 1.03 && n < 57.4) return 1;
    
    return 0;        
}

MODULE = Math::Random::BlumBlumShub	PACKAGE = Math::Random::BlumBlumShub	

PROTOTYPES: DISABLE


void
bbs (outref, p, q, seed, bits_required)
	mpz_t *	outref
	mpz_t *	p
	mpz_t *	q
	mpz_t *	seed
	int	bits_required
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	bbs(outref, p, q, seed, bits_required);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
bbs_seedgen (seed, p, q)
	mpz_t *	seed
	mpz_t *	p
	mpz_t *	q
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	bbs_seedgen(seed, p, q);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

int
monobit (bitstream)
	mpz_t *	bitstream

int
longrun (bitstream)
	mpz_t *	bitstream

int
runs (bitstream)
	mpz_t *	bitstream

int
poker (bitstream)
	mpz_t *	bitstream


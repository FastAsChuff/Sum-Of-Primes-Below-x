#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <assert.h>
#include "/home/simon/Mairsonsprimesieve.c" // https://github.com/FastAsChuff/Primes-List/blob/main/Mairsonsprimesieve.c
#ifndef alignof 
  #define alignof _Alignof
#endif
#define U128 unsigned __int128
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

// See https://math.stackexchange.com/questions/1378286/find-the-sum-of-all-primes-smaller-than-a-big-number/5118835#5118835

// gcc sumofprimes2.c -o sumofprimes2.bin -O3 -Wall -mssse3 -lm

#define NUMOF32BITPRIMES 203280221

// RAM Requirements
// ================
// xi stack 151876932*sizeof(xivals_s) = up to 9.72GB
// Choose cache size...(diminishing returns)
#define SOPXICACHESIZE1 // ~1GB
//#define SOPXICACHESIZE2 // ~3.1GB
//#define SOPXICACHESIZE3 // ~16.6GB
//#define SOPXICACHESIZE4 // ~28.8GB
#ifdef SOPXICACHESIZE1 
#define SOPXICACHETYPE uint32_t
#define SOPXICACHEWIDTH 50000
#define SOPXICACHEHEIGHT 5133
#endif
#ifdef SOPXICACHESIZE2 
#define SOPXICACHETYPE uint32_t
#define SOPXICACHEWIDTH 90000
#define SOPXICACHEHEIGHT 8713
#endif
#ifdef SOPXICACHESIZE3 
#define SOPXICACHETYPE uint64_t
#define SOPXICACHEWIDTH 150000
#define SOPXICACHEHEIGHT 13848
#endif
#ifdef SOPXICACHESIZE4 
#define SOPXICACHETYPE uint64_t
#define SOPXICACHEWIDTH 200000
#define SOPXICACHEHEIGHT 17984
#endif

uint32_t isqrt(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

_Bool inascarrayu32(uint32_t n, uint32_t *array, uint64_t arraysize, uint64_t *index) {
  // Binary Search. Array is assumed sorted ascending. index largest st array[index] <= n.
  if ((n < array[0]) || (n > array[arraysize-1])) {
    *index = -1;
    return false;
  }
  uint64_t leftix = 0;
  uint64_t rightix = arraysize-1;
  if (array[leftix] == n) {
    *index = leftix;
    return true;
  }
  if (array[rightix] == n) {  
    *index = rightix;
    return true;
  }
  uint64_t middleix = (leftix + rightix)/2;
  while (true) {
    if (array[middleix] == n) {
      *index = middleix;
      return true;
    }
    if (array[middleix] < n) {
      leftix = middleix;
    } else {
      rightix = middleix;
    }
    middleix = (leftix + rightix)/2;
    if (middleix == leftix) {
      *index = middleix;
      return array[middleix] == n;
    }
  } 
  assert(true == false); 
  return false;
}

void printU128(U128 n) {
  if (n == 0) {
    printf("0");
    return;
  }
  char buf[40] = {0};
  int bufix = 1;
  while (n) {
    buf[bufix++] = n % 10;
    n /= 10;
  }
  bufix--;
  while (bufix) printf("%c", '0' + buf[bufix--]);
}

typedef struct {
  U128 val1, val2; // xi(x, m) = va1 - p_m*val2 if both are not -1.
  uint64_t x; // x > 0 if element values can be used
  uint32_t m, isqrtx, piisqrtx; // m <= piisqrtx
} xivals_s;

void getxivals(xivals_s *xivals, uint64_t stacksize, uint32_t *primes, uint64_t *primesums, uint32_t numprimes, SOPXICACHETYPE *xicachesmalls) {
  uint64_t xivalsix = 0;
  while (true) {
    _Bool done = false; 
    if (xivals[xivalsix].m == 0) {
      xivals[xivalsix].val1 = ((U128)(xivals[xivalsix].x+1)*xivals[xivalsix].x)/2;
      xivals[xivalsix].val2 = 0;
    }    
    if ((xivals[xivalsix].x <= SOPXICACHEWIDTH) && (xivals[xivalsix].m <= SOPXICACHEHEIGHT)) {
      if (xicachesmalls[xivals[xivalsix].x*SOPXICACHEHEIGHT + xivals[xivalsix].m] != (SOPXICACHETYPE)-1) {
        xivals[xivalsix].val1 = xicachesmalls[xivals[xivalsix].x*SOPXICACHEHEIGHT + xivals[xivalsix].m];
        xivals[xivalsix].val2 = 0;
      }
    } 
    if (xivals[xivalsix].x <= primes[numprimes-1]) {
      if (xivals[xivalsix].m == xivals[xivalsix].piisqrtx) {
        uint64_t pix = 0;
        if (xivals[xivalsix].x >= 2) {
          inascarrayu32(xivals[xivalsix].x, primes, numprimes, &pix);
          pix++;
        }
        xivals[xivalsix].val1 = 1;
        if (pix) xivals[xivalsix].val1 += primesums[pix-1];
        if (xivals[xivalsix].piisqrtx) xivals[xivalsix].val1 -= primesums[xivals[xivalsix].piisqrtx-1];
        xivals[xivalsix].val2 = 0;
      }
    }   
    if (xivals[xivalsix].val1 == (U128)-1) {
      if (xivals[1+xivalsix].x) {
        //assert((xivals[xivalsix].x == xivals[1+xivalsix].x) && (xivals[xivalsix].m == 1+xivals[1+xivalsix].m));
        //assert((xivals[1+xivalsix].val1 != (U128)-1) && (xivals[1+xivalsix].val2 != (U128)-1));
        if (xivals[1+xivalsix].m) {
          xivals[xivalsix].val1 = xivals[1+xivalsix].val1 - primes[xivals[1+xivalsix].m-1]*xivals[1+xivalsix].val2;
        } else {
          xivals[xivalsix].val1 = xivals[1+xivalsix].val1;
        }
        xivals[1+xivalsix].x = 0;
      } else {
        xivals[1+xivalsix].x = xivals[xivalsix].x;
        xivals[1+xivalsix].m = xivals[xivalsix].m-1;
        xivals[1+xivalsix].isqrtx = xivals[xivalsix].isqrtx;
        xivals[1+xivalsix].piisqrtx = xivals[xivalsix].piisqrtx;
        xivalsix++;
        xivals[xivalsix].val1 = (U128)-1; 
        xivals[xivalsix].val2 = (U128)-1; 
        done = true;
      }
    }
    if (!done && (xivals[xivalsix].val2 == (U128)-1) && (xivals[xivalsix].val1 != (U128)-1)) {
      if (xivals[1+xivalsix].x) {
        //assert(xivals[xivalsix].x/primes[xivals[xivalsix].m-1] == xivals[1+xivalsix].x);
        //assert((xivals[1+xivalsix].val1 != (U128)-1) && (xivals[1+xivalsix].val2 != (U128)-1));
        if (xivals[1+xivalsix].m) {
          xivals[xivalsix].val2 = xivals[1+xivalsix].val1 - primes[xivals[1+xivalsix].m-1]*xivals[1+xivalsix].val2;
        } else {
          xivals[xivalsix].val2 = xivals[1+xivalsix].val1;
        }
    // s(x) = xi(x, pi(a)) + s(a) - 1 if a >= isqrt(x)
    // s(x) = xi(x, m) + s(p_m) - 1 if p_m >= isqrt(x)
    // s(x) = xi(x, m) + s(p_m) - 1 if m >= piisqrtx
    // s(x) = xi(x, m) + s(p_m) - 1 = xi(x, piisqrtx) + s(isqrtx) - 1 if m >= piisqrtx
    // xi(x, m) = xi(x, piisqrtx) + s(isqrtx) - s(p_m) if m >= piisqrtx
        if (1+xivals[1+xivalsix].m < xivals[xivalsix].m) {
          //assert(xivals[1+xivalsix].m == xivals[1+xivalsix].piisqrtx);
          if (xivals[1+xivalsix].piisqrtx) {
            xivals[xivalsix].val2 = (xivals[xivalsix].val2 + primesums[xivals[1+xivalsix].piisqrtx-1]) - primesums[xivals[xivalsix].m-2];
          } else {
            xivals[xivalsix].val2 -= primesums[xivals[xivalsix].m-2];
          }
        }
        xivals[1+xivalsix].x = 0;
      } else {
        xivals[1+xivalsix].x = xivals[xivalsix].x/primes[xivals[xivalsix].m-1];
        if (xivals[1+xivalsix].x) {
          xivals[1+xivalsix].isqrtx = isqrt(xivals[1+xivalsix].x);
          uint64_t piisqrtx = 0;
          if (xivals[1+xivalsix].isqrtx >= 2) {
            inascarrayu32(xivals[1+xivalsix].isqrtx, primes, numprimes, &piisqrtx);
            piisqrtx++;
          }
          xivals[1+xivalsix].piisqrtx = piisqrtx;
          xivals[1+xivalsix].m = MIN(xivals[xivalsix].m-1, piisqrtx);
          xivalsix++;
          xivals[xivalsix].val1 = (U128)-1; 
          xivals[xivalsix].val2 = (U128)-1; 
          done = true;
        } else {
          xivals[xivalsix].val2 = 0;
        }
      }
    }
    /*
    printf("%lu: x=%lu m=%u", xivalsix, xivals[xivalsix].x, xivals[xivalsix].m);
    if (xivals[xivalsix].val1 != (U128)-1) {
      printf(" val1=");
      printU128(xivals[xivalsix].val1);
    }
    if (xivals[xivalsix].val2 != (U128)-1) {
      printf(" val2=");
      printU128(xivals[xivalsix].val2);
    }
    if ((xivals[xivalsix].val1 != (U128)-1) && (xivals[xivalsix].val2 != (U128)-1)) {
      printf(" val=");
      U128 val = xivals[xivalsix].val1;
      if (xivals[xivalsix].m) val -= primes[xivals[xivalsix].m-1]*xivals[xivalsix].val2;
      printU128(val);
    }
    printf("\n");
    sleep(1);
    */
    if ((xivals[xivalsix].val2 != (U128)-1) && (xivals[xivalsix].val1 != (U128)-1)) {
      if (xivalsix == 0) break;
      if ((xivals[xivalsix].x <= SOPXICACHEWIDTH) && (xivals[xivalsix].m <= SOPXICACHEHEIGHT)) {
        if (xicachesmalls[xivals[xivalsix].x*SOPXICACHEHEIGHT + xivals[xivalsix].m] == (SOPXICACHETYPE)-1) {
          U128 val = xivals[xivalsix].val1;
          if (xivals[xivalsix].m) val -= primes[xivals[xivalsix].m-1]*xivals[xivalsix].val2;
          xicachesmalls[xivals[xivalsix].x*SOPXICACHEHEIGHT + xivals[xivalsix].m] = val;
        }  
      }
      xivalsix--;      
      //assert(((xivals[xivalsix].val1 != (U128)-1) && (xivals[xivalsix].val2 == (U128)-1)) || (xivals[xivalsix].val1 == (U128)-1));
    }
  }
}

U128 getxi(uint64_t x, uint32_t m, uint32_t **primes, uint64_t **primesums, uint32_t *numprimes, SOPXICACHETYPE *xicachesmalls) {
  if (m == 0) return ((U128)(x+1)*x)/2;
  if (x == 0) return 0;
  if (m >= x) return 1;
  if (m <= *numprimes) {
    if ((*primes)[m-1] >= x) return 1;
  }
  if ((x <= SOPXICACHEWIDTH) && (m <= SOPXICACHEHEIGHT)) {
    if (xicachesmalls[x*SOPXICACHEHEIGHT + m] != (SOPXICACHETYPE)-1) return xicachesmalls[x*SOPXICACHEHEIGHT + m];
  }
  uint32_t isqrtx = isqrt(x);
  uint64_t piisqrtx = 0;
  if (isqrtx >= 2) {
    inascarrayu32(isqrtx, *primes, *numprimes, &piisqrtx);
    piisqrtx++;
  }
  if (piisqrtx) {
    if (m > piisqrtx) return (getxi(x, piisqrtx, primes, primesums, numprimes, xicachesmalls) + (*primesums)[piisqrtx-1]) - (*primesums)[m-1];
  }
  xivals_s *xivals = aligned_alloc(alignof(xivals_s), (m+7)*sizeof(xivals_s));
  assert(xivals); 
  for(uint64_t i=0; i<m+7; i++) xivals[i].x = 0;
  xivals[0].val1 = (U128)-1; 
  xivals[0].val2 = (U128)-1; 
  xivals[0].x = x; 
  xivals[0].m = m; 
  xivals[0].isqrtx = isqrtx; 
  xivals[0].piisqrtx = piisqrtx;
  getxivals(xivals, m, *primes, *primesums, *numprimes, xicachesmalls);  
  U128 xi = xivals[0].val1 - (*primes)[m-1]*xivals[0].val2;
  free(xivals);
  return xi;
}

uint64_t *makeprimesums(uint32_t *primes, uint32_t numprimes) {
  uint64_t sum = 0;
  uint64_t *primesums = malloc(numprimes*(sizeof(uint64_t)));
  if (primesums) {
    for (uint32_t i=0; i<numprimes; i++) primesums[i] = (sum += primes[i]);
  }
  return primesums;
}

U128 sumofprimes(uint64_t x, uint32_t **primes, uint64_t **primesums, uint32_t *numprimes, SOPXICACHETYPE *xicachesmalls) {
  if (x < 2) return 0;
  if (x <= (*primes)[*numprimes-1]) {
    uint64_t pix;
    inascarrayu32(x, *primes, *numprimes, &pix);
    return (*primesums)[pix];
  }
  uint32_t isqrtx = isqrt(x);
  if ((isqrtx > (*primes)[*numprimes - 1]) || ((isqrtx > 5000000U) && (*numprimes < NUMOF32BITPRIMES))) {
    free(*primes);
    free(*primesums);
    if (isqrtx > 5000000U) {
      *primes = Mairsonsprimesieve(0xffffffffU, numprimes);
    } else {
      *primes = Mairsonsprimesieve(2*isqrtx, numprimes);
    }
    assert(*primes);
    *primesums = makeprimesums(*primes, *numprimes);
    assert(*primesums);
  }
  uint64_t piisqrtx;
  inascarrayu32(isqrtx, *primes, *numprimes, &piisqrtx);
  piisqrtx++;
  //printf("There are %lu primes <= %u.\n", piisqrtx, isqrtx);
  //exit(0);
  return getxi(x, piisqrtx, primes, primesums, numprimes, xicachesmalls) + sumofprimes(isqrtx, primes, primesums, numprimes, xicachesmalls) - 1;
}

int main(int argc, char* argv[]) {
  uint64_t maxMB = 10;
  if (argc > 1) {
    maxMB = atol(argv[1]);
    maxMB = (maxMB >= 10 ? maxMB : 10);
  }
  uint32_t numprimes;
  uint32_t *primes = Mairsonsprimesieve(100000000u, &numprimes);
  assert(primes);
  uint64_t *primesums = makeprimesums(primes, numprimes);
  uint64_t xicachesmallsbytes = (1+SOPXICACHEWIDTH)*((1+SOPXICACHEHEIGHT)*sizeof(SOPXICACHETYPE));
  SOPXICACHETYPE *xicachesmalls = aligned_alloc(alignof(SOPXICACHETYPE), xicachesmallsbytes);
  assert(xicachesmalls);
  memset(xicachesmalls, 0xffu, xicachesmallsbytes);
  printf("This program calculates the sum of all the prime numbers less than or equal to input value x <= 10^19.\nAuthor: Simon Goater Jan 2026\nRAM requirement approx. %f + %lu*pi(sqrt(x))/1000000000 GB.\n", 
(((double)sizeof(SOPXICACHETYPE)*SOPXICACHEWIDTH)*SOPXICACHEHEIGHT + (double)NUMOF32BITPRIMES*(sizeof(uint32_t) + sizeof(uint64_t)))/1000000000.0, sizeof(xivals_s));
  while (true) {
    uint64_t x;
    printf("Enter x: ");
    fflush(stdout);
    assert(1 == scanf("%lu", &x));
    assert(x <= 10000000000000000000ULL);
    printf("Calculating sum of primes less than or equal to %lu.\n", x);  
    uint64_t timestart = time(0);
    U128 result = sumofprimes(x, &primes, &primesums, &numprimes, xicachesmalls);
    uint64_t timeend = time(0);
    printf("The sum of primes less than or equal to %lu is ", x); 
    printU128(result);
    printf(".\n");
    uint64_t timesecs = timeend - timestart;
    if (timesecs > 5) printf("%f Gx/s\n", x/(1000000000.0*timesecs));
    uint64_t timehrs = timesecs / 3600;
    uint64_t timemins = (timesecs - (3600*timehrs)) / 60;
    timesecs -= 60*(timemins + 60*timehrs);
    printf("Time: ");
    if (timehrs) printf("%lu hrs ", timehrs);
    if (timehrs || timemins) printf("%lu mins ", timemins);
    printf("%lu secs\n", timesecs);
  }
}

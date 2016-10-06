#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <mm_malloc.h>
#include <x86intrin.h>

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) < (B)) ? (B) : (A))

#define L1D_CACHE_SIZE 32768
#define ALIGN_FACTOR 64

typedef unsigned int u32;
typedef unsigned long long u64;
typedef unsigned char u8;

int main(int argc, char* argv[]) {
  struct timespec start, stop;
  u32 count;
  u64 number = 1000000000;
  if (argc > 1)
    number = strtoull(argv[1], NULL, 10);

  clock_gettime(CLOCK_REALTIME, &start);
  const u32 SqrtNum = (u32)ceil(sqrt(number));
  u8 *restrict isComposite = (u8*)_mm_malloc(sizeof(u8) * SqrtNum, ALIGN_FACTOR);
  memset (isComposite, 0, SqrtNum);
  for (u32 i = 2; i * i <= SqrtNum; ++i) {
    if (!isComposite[i])
      for (u32 j = i * i; j <= SqrtNum; j += i)
        isComposite[j] = 1;
  }
  const u32 SegmentSize = MAX(SqrtNum, L1D_CACHE_SIZE);
  u8 *restrict sieve = (u8*)_mm_malloc(sizeof(u8) * SegmentSize, ALIGN_FACTOR);
  u32 *restrict primes = (u32*)_mm_malloc(sizeof(u32) * SqrtNum / 4, ALIGN_FACTOR);
  u32 *restrict nextIndex = (u32*)_mm_malloc(sizeof(u32) * SqrtNum / 4, ALIGN_FACTOR);

  u64 currNum = 3;
  u64 numIndex = 3;
  u32 primeCounter = 0;
  count = (number > 1);

  for (u64 low = 0; low <= number; low += SegmentSize) {
    u64 high = MIN(low + SegmentSize - 1, number);
    memset(sieve, 1, SegmentSize);
    for (; currNum * currNum <= high; currNum += 2)
      if (!isComposite[currNum]) {
        primes[primeCounter] = currNum;
        nextIndex[primeCounter] = currNum * currNum - low;
        ++primeCounter;
      }
    for (u32 i = 0; i < primeCounter; ++i) {
      u32 ni = nextIndex[i];
      const u32 incr = primes[i] << 1; // skip evens
      for (; ni < SegmentSize; ni += incr)
        sieve[ni] = 0; // mark composites
      nextIndex[i] = ni - SegmentSize;
    }
    for (; numIndex + 8 <= high; numIndex += 8)
      count += _mm_popcnt_u64(*((u64*)(sieve + numIndex - low))) - 4;
    for (; numIndex + 4 <= high; numIndex += 4)
      count += _mm_popcnt_u32(*((u32*)(sieve + numIndex - low))) - 2;
    for (; numIndex <= high; numIndex += 2)
      count += sieve[numIndex - low];
  }
  _mm_free (nextIndex);
  _mm_free (primes);
  _mm_free (sieve);
  _mm_free (isComposite);
  clock_gettime(CLOCK_REALTIME, &stop);

  printf("There are %u prime numbers under %llu.\n", count, number);
  printf("Time (ns): %llu\n", 1000000000LL * (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec));
  return 0;
}

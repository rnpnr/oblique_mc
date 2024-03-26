#include <stdint.h>
#include <time.h>

typedef uint64_t u64;

#include "utils.h"

static u64 rand_state[2];

void
random_init(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	rand_state[0] = (intptr_t)&printf ^ ts.tv_sec
	              ^ ((u64)ts.tv_nsec * 0xAC5533CD);
	rand_state[1] = (intptr_t)&malloc ^ ts.tv_sec
	              ^ ((u64)ts.tv_nsec * 0xAC5533CD);
}

static u64
xoroshiro128plus(u64 s[2])
{
	u64 s0 = s[0];
	u64 s1 = s[1];
	u64 result = s0 + s1;
	s1 ^= s0;
	s[0] = ((s0 << 24) | (s0 >> 40)) ^ s1 ^ (s1 << 16);
	s[1] = (s1 << 37) | (s1 >> 27);
	return result;
}

double
random_uniform(void)
{
	return xoroshiro128plus(rand_state) / (double)UINT64_MAX;
}

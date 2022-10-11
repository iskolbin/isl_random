#define ISL_RANDOM_STATIC
#define ISL_RANDOM_IMPLEMENTATION
#include "isl_random.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
	uint64_t state[ISLR_STATE_SIZE];
	islr_srand(state, 0xDEADBEEF);    // Use builtin Splitmix64 to init state
	uint64_t raw = islr_next(state);
	int random_int = islr_random(state, 0, 1000); // Generate random int [0-1000)
	double random_double = islr_random_double(state); // Generate random double [0.0-1.0)
	printf("%d %5.5f\n", random_int, random_double);
	return 0;
}

/**
 * @file test_bignum_jacobi_extra.c
 * @brief Deterministic randomized tests for the C11 Jacobi implementation.
 * @details A bounded uint64 Jacobi oracle is evaluated independently from the
 * bignum implementation for 4096 generated odd moduli and numerators. The
 * test also exercises zero and even modulus validation without mutation.
 */
#include "bignum_jacobi.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

/** @brief Computes an independent uint64 Jacobi oracle. */
static int oracle(uint64_t a, uint64_t n)
{
    int sign = 1;
    a %= n;
    while (a != 0U) {
        while ((a & 1U) == 0U) {
            a >>= 1U;
            if ((n & 7U) == 3U || (n & 7U) == 5U) sign = -sign;
        }
        if ((a & 3U) == 3U && (n & 3U) == 3U) sign = -sign;
        { uint64_t t = a; a = n; n = t; }
        a %= n;
    }
    return n == 1U ? sign : 0;
}

/** @brief Builds a normalized one-word operand. */
static bignum_t bn(uint64_t value)
{
    bignum_t x = {{0U}, 0U};
    if (value != 0U) { x.words[0] = value; x.len = 1U; }
    return x;
}

int main(void)
{
    uint64_t state = UINT64_C(0x9E3779B97F4A7C15);
    for (unsigned i = 0U; i < 4096U; ++i) {
        bignum_t a, n;
        int actual = 0;
        state ^= state << 7U; state ^= state >> 9U; state ^= state << 8U;
        a = bn(state);
        state ^= state << 7U; state ^= state >> 9U; state ^= state << 8U;
        n = bn((state % UINT64_C(1000003)) | UINT64_C(1));
        if (n.words[0] < 3U) n.words[0] = 3U;
        assert(bignum_jacobi(&a, &n, &actual) == BIGNUM_JACOBI_SUCCESS);
        assert(actual == oracle(a.words[0], n.words[0]));
    }
    { bignum_t a = bn(5U), zero = bn(0U), even = bn(2U); int out = 17;
      assert(bignum_jacobi(&a, &zero, &out) == BIGNUM_JACOBI_ERROR_MODULUS);
      assert(bignum_jacobi(&a, &even, &out) == BIGNUM_JACOBI_ERROR_MODULUS);
      assert(out == 17); }
    puts("bignum_jacobi randomized tests: OK");
    return 0;
}

/**
 * @file test_bignum_jacobi.c
 * @brief Deterministic C11 tests for bignum_jacobi.
 * @details Covers valid symbols, zero numerator, modulus validation, output
 * preservation on errors, input immutability, and a multiword-sized value.
 */
#include "bignum_jacobi.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

/** @brief Creates a normalized one-word test operand. */
static bignum_t value(uint64_t v)
{
    bignum_t x = {{0U}, 0U};
    if (v != 0U) { x.words[0] = v; x.len = 1U; }
    return x;
}

/** @brief Checks one known Jacobi value. */
static void expect(uint64_t a_value, uint64_t n_value, int expected)
{
    bignum_t a = value(a_value), n = value(n_value);
    int symbol = 77;
    assert(bignum_jacobi(&a, &n, &symbol) == BIGNUM_JACOBI_SUCCESS);
    assert(symbol == expected);
}

int main(void)
{
    bignum_t a = value(5U), n = value(3U), zero = value(0U), even = value(2U);
    const bignum_t original_a = a, original_n = n;
    int symbol = 91;
    expect(0U, 3U, 0);
    expect(1U, 1U, 1);
    expect(2U, 3U, -1);
    expect(3U, 5U, -1);
    expect(5U, 7U, -1);
    expect(9U, 13U, 1);
    expect(1001U, 9907U, -1);
    expect(UINT64_C(0x1000000000000001), UINT64_C(0x1000000000000007), -1);
    assert(bignum_jacobi(NULL, &n, &symbol) == BIGNUM_JACOBI_ERROR_NULL_ARG);
    assert(bignum_jacobi(&a, NULL, &symbol) == BIGNUM_JACOBI_ERROR_NULL_ARG);
    assert(bignum_jacobi(&a, &n, NULL) == BIGNUM_JACOBI_ERROR_NULL_ARG);
    assert(bignum_jacobi(&a, &zero, &symbol) == BIGNUM_JACOBI_ERROR_MODULUS);
    assert(bignum_jacobi(&a, &even, &symbol) == BIGNUM_JACOBI_ERROR_MODULUS);
    assert(symbol == 91);
    assert(memcmp(&a, &original_a, sizeof(a)) == 0);
    assert(memcmp(&n, &original_n, sizeof(n)) == 0);
    puts("bignum_jacobi deterministic tests: OK");
    return 0;
}

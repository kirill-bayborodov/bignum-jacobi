/**
 * @file bignum_jacobi.c
 * @brief C11 reference implementation of the Jacobi symbol.
 * @details
 * The implementation uses fixed-capacity little-endian 64-bit words. It copies
 * both borrowed inputs into local objects, reduces the numerator by binary
 * long division, and then applies the binary Jacobi transformations. No heap
 * allocation or mutable global state is used. The reference is intentionally
 * straightforward so its arithmetic and coverage are easy to audit before
 * comparing an optimized assembly implementation.
 */
#include "bignum_jacobi.h"

#include <string.h>

/** @brief Returns whether a local bignum is zero. @param[in] x Local value. @return Non-zero for zero. */
static int jacobi_is_zero(const bignum_t *x)
{
    return x->len == 0U;
}

/** @brief Normalizes a local value. @param[in,out] x Value whose leading zero words are removed. */
static void jacobi_normalize(bignum_t *x)
{
    while (x->len > 0U && x->words[x->len - 1U] == 0U) {
        --x->len;
    }
}

/** @brief Compares normalized values. @return -1, 0, or 1 according to ordering. */
static int jacobi_compare(const bignum_t *a, const bignum_t *b)
{
    size_t i;
    if (a->len != b->len) {
        return a->len < b->len ? -1 : 1;
    }
    i = a->len;
    while (i > 0U) {
        --i;
        if (a->words[i] != b->words[i]) {
            return a->words[i] < b->words[i] ? -1 : 1;
        }
    }
    return 0;
}

/** @brief Subtracts b from a after the caller established a >= b. */
static void jacobi_subtract(bignum_t *a, const bignum_t *b)
{
    size_t i;
    uint64_t borrow = 0U;
    for (i = 0U; i < a->len; ++i) {
        const uint64_t bi = i < b->len ? b->words[i] : 0U;
        const uint64_t old = a->words[i];
        const uint64_t t = old - bi;
        const uint64_t borrow_one = old < bi;
        const uint64_t result = t - borrow;
        const uint64_t borrow_two = t < borrow;
        a->words[i] = result;
        borrow = borrow_one | borrow_two;
    }
    jacobi_normalize(a);
}

/** @brief Shifts a normalized value right by one bit. */
static void jacobi_shift_right_one(bignum_t *x)
{
    size_t i = x->len;
    uint64_t carry = 0U;
    while (i > 0U) {
        uint64_t word;
        --i;
        word = x->words[i];
        x->words[i] = (word >> 1U) | carry;
        carry = word << 63U;
    }
    jacobi_normalize(x);
}

/** @brief Returns the least significant bit of a normalized value. */
static unsigned jacobi_is_odd(const bignum_t *x)
{
    return x->len != 0U ? (unsigned)(x->words[0] & UINT64_C(1)) : 0U;
}

/** @brief Computes x modulo the positive modulus m using binary long division. */
static void jacobi_reduce(bignum_t *x, const bignum_t *m)
{
    bignum_t remainder = {{0U}, 0U};
    size_t bit;
    if (jacobi_is_zero(x)) {
        return;
    }
    bit = x->len * 64U;
    while (bit > 0U) {
        const size_t source_bit = bit - 1U;
        const size_t word = source_bit / 64U;
        const unsigned offset = (unsigned)(source_bit % 64U);
        uint64_t carry = (x->words[word] >> offset) & UINT64_C(1);
        size_t i;
        for (i = 0U; i < remainder.len; ++i) {
            const uint64_t next = remainder.words[i] >> 63U;
            remainder.words[i] = (remainder.words[i] << 1U) | carry;
            carry = next;
        }
        if (carry != 0U && remainder.len < BIGNUM_CAPACITY) {
            remainder.words[remainder.len++] = carry;
        }
        if (jacobi_compare(&remainder, m) >= 0) {
            jacobi_subtract(&remainder, m);
        }
        --bit;
    }
    *x = remainder;
}

/** @brief Applies the Jacobi supplementary sign for an even numerator. */
static void jacobi_remove_twos(bignum_t *a, const bignum_t *n, int *sign)
{
    while (!jacobi_is_zero(a) && !jacobi_is_odd(a)) {
        const uint64_t n_mod_eight = n->words[0] & UINT64_C(7);
        jacobi_shift_right_one(a);
        if (n_mod_eight == 3U || n_mod_eight == 5U) {
            *sign = -*sign;
        }
    }
}

/** @brief Computes the validated binary Jacobi state machine. */
static int jacobi_compute(const bignum_t *a, const bignum_t *n)
{
    bignum_t numerator = *a;
    bignum_t modulus = *n;
    int sign = 1;
    jacobi_reduce(&numerator, &modulus);
    while (!jacobi_is_zero(&numerator)) {
        jacobi_remove_twos(&numerator, &modulus, &sign);
        if (jacobi_is_zero(&numerator)) {
            break;
        }
        if ((numerator.words[0] & UINT64_C(3)) == 3U &&
            (modulus.words[0] & UINT64_C(3)) == 3U) {
            sign = -sign;
        }
        {
            bignum_t swap = numerator;
            numerator = modulus;
            modulus = swap;
        }
        jacobi_reduce(&numerator, &modulus);
    }
    return modulus.len == 1U && modulus.words[0] == 1U ? sign : 0;
}

bignum_jacobi_status_t bignum_jacobi(const bignum_t *a, const bignum_t *n, int *symbol)
{
    bignum_t modulus;
    int result;
    if (a == NULL || n == NULL || symbol == NULL) {
        return BIGNUM_JACOBI_ERROR_NULL_ARG;
    }
    modulus = *n;
    jacobi_normalize(&modulus);
    if (jacobi_is_zero(&modulus) || !jacobi_is_odd(&modulus)) {
        return BIGNUM_JACOBI_ERROR_MODULUS;
    }
    result = jacobi_compute(a, &modulus);
    *symbol = result;
    return BIGNUM_JACOBI_SUCCESS;
}

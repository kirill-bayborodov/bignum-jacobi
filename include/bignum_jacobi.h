/**
 * @file bignum_jacobi.h
 * @brief Public API for the Jacobi symbol of bounded unsigned integers.
 * @details
 * The module evaluates (a/n) for non-negative bignum_t operands using the
 * standard binary Jacobi algorithm. The modulus must be positive and odd.
 * Inputs are caller-owned borrowed values and are never modified. The result
 * is written to caller-owned storage only after successful validation and is
 * one of -1, 0, or 1. The operation performs no allocation and is safe for
 * concurrent calls when callers do not concurrently mutate the input objects.
 * The optimized implementation preserves the System V AMD64 ABI and has the
 * same status and output contract as the C11 reference implementation.
 * @warning The function does not interpret a zero or even modulus as a valid
 * Jacobi modulus; such input returns a named error and leaves the output
 * unchanged. The input bignum_t values must be normalized.
 */
#ifndef BIGNUM_JACOBI_H
#define BIGNUM_JACOBI_H

#include <bignum.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Reports validation and evaluation status for bignum_jacobi.
 * @details A successful status guarantees that the output symbol was written.
 * Every error status leaves the caller-owned output object unchanged and is
 * suitable for correction and retry with valid operands.
 */
typedef enum bignum_jacobi_status {
    BIGNUM_JACOBI_SUCCESS = 0, /**< Symbol computed; output is exactly -1, 0, or 1. */
    BIGNUM_JACOBI_ERROR_NULL_ARG = -1, /**< An input or output pointer is NULL; output is unchanged. */
    BIGNUM_JACOBI_ERROR_MODULUS = -2 /**< Modulus is zero or even; output is unchanged and no symbol is defined. */
} bignum_jacobi_status_t;

/**
 * @brief Computes the Jacobi symbol of two non-negative bounded integers.
 * @details The binary algorithm repeatedly removes factors of two from the
 * numerator, applies the supplementary law for 2, swaps numerator/modulus,
 * applies quadratic reciprocity, and reduces the numerator modulo the modulus.
 * The result is 0 when the final numerator is zero and the remaining modulus
 * is greater than one; otherwise it is the accumulated sign when the modulus
 * reaches one. Time complexity is O(W^2) in the fixed capacity W-word model;
 * temporary storage is O(W) and is stack-local.
 * @param[in] a Borrowed normalized non-negative numerator; it is not modified.
 * @param[in] n Borrowed normalized positive odd modulus; it is not modified.
 * @param[out] symbol Caller-allocated integer receiving -1, 0, or 1 on success.
 *                    It may alias neither input and remains unchanged on error.
 * @return BIGNUM_JACOBI_SUCCESS, BIGNUM_JACOBI_ERROR_NULL_ARG, or
 *         BIGNUM_JACOBI_ERROR_MODULUS.
 * @pre `a`, `n`, and `symbol` are valid pointers; `a` and `n` are normalized.
 * @post Inputs are unchanged. On success `*symbol` is a Jacobi value.
 * @warning The modulus must be odd; this API does not silently normalize it.
 * @thread_safety Safe for independent immutable inputs and distinct outputs.
 */
bignum_jacobi_status_t bignum_jacobi(const bignum_t *a, const bignum_t *n, int *symbol);

#ifdef __cplusplus
}
#endif

#endif /* BIGNUM_JACOBI_H */

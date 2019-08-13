#ifndef MATHUTIL_H
#define MATHUTIL_H

/**
 * @file mathutil.h
 * @author Paulo Fonseca
 * 
 * @brief Math utility functions.
 */


/**
 * @brief Computes the generalised multiplicative floor function defined as
 * the integral multiplicative factor corresponding to largest multiple of a
 * given fixed (finite) base no greater than a given (finite) threshold value.
 * If the base is zero, the function is ill-defined since all its multiples
 * will be equally zero and hence either all of them will be no greater
 * than the threshold, or none of them will. In the latter case, the function
 * shall be define as -INFINITY whereas in the former case it will be defined
 * as zero. In the other cases, the function is well-defined and has
 * a finite value (possibly negative, if value and base have opposite signs).
 * Note that mult_floor(n, 1)==floor(n).
 *
 * @param value  The limiting finite value.
 * @param base The fixed finite base.
 * @return  The integral multiplicative factor corresponding to largest
 * multiple of "base" no greater than "value".
 */
double multfloor(double value, double base);

/**
 * @brief Computes the generalised multiplicative ceil function defined as
 * the integral multiplicative factor corresponding to smallest multiple of a
 * given fixed (finite) base no smaller than a given (finite) threshold value.
 * If the base is zero, the function is ill-defined since all its multiples
 * will be equally zero and hence either all of them will be no smaller
 * than the threshold, or none of them will. In the latter case, the function
 * shall be define as +INFINITY whereas in the former case it will be defined
 * as zero. In the other cases, the function is well-defined and has
 * a finite value (possibly negative, if value and base have opposite signs).
 * Note that mult_ceil(n, 1)==ceil(n).
 *
 * @param value  The limiting finite value.
 * @param base The fixed finite base.
 * @return  The integral multiplicative factor corresponding to smallest
 * multiple of "base" no smaller than "value".
 */
double multceil(double val, double fact);


#endif
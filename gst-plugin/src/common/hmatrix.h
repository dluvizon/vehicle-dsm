/** @file hmatrix.h
 * @author Diogo Luvizon <diogo@luvizon.com
 * @date 26/08/2014
 */

#ifndef __HMATRIX_H
#define __HMATRIX_H

/** @struct hmatrix
 * @brief Defines a simple 3x3 homography matrix
 *
 *  a |--|--|
 *  b |--|--|
 *  c |--|--|
 *    0  1  2
 */
struct hmatrix {
	float a0, a1, a2;
	float b0, b1, b2;
	float c0, c1, c2;
};

#endif /* __HMATRIX_H */

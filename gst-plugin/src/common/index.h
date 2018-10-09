/** Lightweight Computer Vision Library
 * ------------------------------------
 * @file lwcv-index.h
 * @brief Header file for index operations in 2D images.
 * @copyright <TODO>
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @version 0.1
 * @date 23/06/2014
 * @section Detailed Description
 * <pre>
 * Definitions for the border pixels around Mi in a frame
 *
 *                            R O W
 * 	------    --------------------------
 * 	| -2 | .. | -2 | -1 | +0 | +1 | +2 | .. 
 * 	------    --------------------------
 * 	| -1 | .. | -2 | -1 | +0 | +1 | +2 | ..
 * 	------    ----------======----------     L
 * 	| +0 | .. | -2 | -1 # Mi # +1 | +2 | ..  I
 * 	------    ----------======----------     N
 * 	| +1 | .. | -2 | -1 | +0 | +1 | +2 | ..  E
 * 	------    -------------------------- 
 * 	| +2 | .. | -2 | -1 | +0 | +1 | +2 | ..
 * 	------    --------------------------
 *
 * Use the macro INDEX to handler the index.
 * </pre>
 */

#ifndef __LWCV_INDEX_H
#define __LWCV_INDEX_H

/** This macro handle the pixel index around a master index - mi.
 * @param mi Master index.
 * @param width Image width.
 * @param row The row number relative to the master index.
 * @param line The line number relative to the master index.
 * @see Detailed Description
 */
#define INDEX(mi, width, row, line) \
	((int)(mi) + (int)(row) + ((int)(line) * (int)(width)))

#define INDEX_INIT(width, row, line) \
	INDEX(0, width, row, line)

#endif /* __LWCV_INDEX_H */

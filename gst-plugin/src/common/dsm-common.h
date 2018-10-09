/** @file dsm-common.h
 * @author Diogo Luvizon <diogo@luvizon.com
 * @date 20/10/2014
 */

#ifndef __DSM_COMMON_H
#define __DSM_COMMON_H

#include <stdint.h>

#define NUM_OF_STILL_BUF 2

/** @struct dsm_common_data
 * @brief Defines the common data (shared) between the gstreamer plugin part
 * and the vehicle-dsm module.
 */
struct dsm_common_data {
	const char *output_path;	/**< Path to output folder. */
	const char *groundtruh_file;	/**< Path to groundtruth file. */
	const char *vehicles_xml;	/**< Path to vehicles.xml file. */
	uint8_t *vdata[NUM_OF_STILL_BUF];	/**< Vector of buffers. */
	int icurr;		/**< Index of the current buffer. */
	int nhold;		/**< Number of hold buffers. */
	long int iframe;	/**< Index of the current frame from video. */

	/* FIXME: verify when the timestamp will overflow. */
	double timestamp;
};

/**
 * Returns a pointer to the past buffer indexed by index.
 * @param cd Pointer for a dsm_common_data structure filled with buffers.
 * @param index Index of the desired buffer, e.g. 0 for the current buffer,
 * 1 for the first past buffer, 2, ..., NUM_OF_STILL_BUF - 1. If the index is
 * greater than NUM_OF_STILL_BUF - 1, it will return NULL.
 * @return Pointer to the desired buffer.
 */
inline uint8_t *grab_past_buffer(struct dsm_common_data *cd, int index)
{
	if (index >= NUM_OF_STILL_BUF) {
		return NULL;
	}
	return cd->vdata[(index + cd->icurr) % NUM_OF_STILL_BUF];
}

#endif /* __DSM_COMMON_H */

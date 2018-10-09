/**
 * @file vehicles.h
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 18/11/2014
 */

#ifndef __VEHICLES_H
#define __VEHICLES_H

#include <stdbool.h>
#include <stdint.h>

#include "xml/pugixml.hpp"
#include "config.h"

struct vehicle_mark {
	int iframe;
	int lane;
	bool plate;
	bool sema;
	bool moto;
	bool radar;
	bool used;
	int x;
	int y;
	int w;
	int h;
	int frame_start;
	int frame_end;
	double speed;
};

struct vehicle_table {
	struct vehicle_mark **vehicles;
	size_t size;
};

/**
 * Create a new vehicle_table structure and initialize with the given file.
 * @param vehicles_xml Path to a xml file with the plates of the vehicles.
 * @return A newly and filled vehicle_table structure.
 */
struct vehicle_table *gtruth_vehicles_new(const char *vehicles_xml);

void gtruth_print_vehicles(const struct vehicle_table *vt);

#endif /* __VEHICLES_H */

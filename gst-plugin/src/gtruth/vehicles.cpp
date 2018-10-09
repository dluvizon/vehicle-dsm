/**
 * @file vehicles.cpp
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 18/11/2014
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#include "vehicles.h"
#include "font.h"

#define _DLEVEL GROUNDTRUTH_DLEVEL
#include "debug.h"

static void parse_vehicles_xml(struct vehicle_table *vt, pugi::xml_node &gt)
{
	pugi::xml_node vehicle;
	pugi::xml_node region;
	pugi::xml_node radar;
	struct vehicle_mark *v;
	int iframe;

	for (vehicle = gt.first_child(); vehicle;
			vehicle = vehicle.next_sibling()) {
		iframe = atoi(vehicle.attribute("iframe").value());
		if (iframe >= vt->size) {
			/* something is wrong, goes out */
			break;
		}
		v = (struct vehicle_mark *) malloc(sizeof(struct vehicle_mark));
		v->used = false;
		v->iframe = iframe;
		v->lane = atoi(vehicle.attribute("lane").value());
		if (strstr(vehicle.attribute("plate").value(), "True")) {
			v->plate = true;
		} else {
			v->plate = false;
		}
		if (strstr(vehicle.attribute("sema").value(), "True")) {
			v->sema = true;
		} else {
			v->sema = false;
		}
		if (strstr(vehicle.attribute("moto").value(), "True")) {
			v->moto = true;
		} else {
			v->moto = false;
		}
		if (strstr(vehicle.attribute("radar").value(), "True")) {
			v->radar = true;
		} else {
			v->radar = false;
		}
		region = vehicle.child("region");
		v->x = atoi(region.attribute("x").value());
		v->y = atoi(region.attribute("y").value());
		v->w = atoi(region.attribute("w").value());
		v->h = atoi(region.attribute("h").value());
		if (v->radar) {
			radar = vehicle.child("radar");
			v->frame_start =
				atoi(radar.attribute("frame_start").value());
			v->frame_end =
				atoi(radar.attribute("frame_end").value());
			v->speed = atof(radar.attribute("speed").value());
		}
		/* store this vehicle in vehicle_table at position iframe */
		if (0 != vt->vehicles[iframe]) {
			print_warn("Damn, the index %d in vehicle_table "
					"is not empty, now skipping", iframe);
		} else {
			vt->vehicles[iframe] = v;
		}
	}
}

struct vehicle_table *gtruth_vehicles_new(const char *vehicles_xml)
{
	/* for handle XML document */
	pugi::xml_document doc;
	pugi::xml_parse_result result;
	pugi::xml_node root;
	pugi::xml_node gtruth;
	pugi::xml_node videoframes;
	struct vehicle_table *vt;
	size_t nbytes;
	int i;

	/* create a new vehicle_table structure */
	vt = new struct vehicle_table;

	/* initialize the xml lib with the given file */
	result = doc.load_file(vehicles_xml);
	root = doc.child("GroundTruthRoot");
	gtruth = root.child("gtruth");
	videoframes = root.child("videoframes");

	vt->size = atoi(videoframes.attribute("total").value());
	nbytes = vt->size * sizeof(struct vehicle_mark *);
	vt->vehicles = (struct vehicle_mark **) malloc(nbytes);
	memset(vt->vehicles, 0, nbytes);
	parse_vehicles_xml(vt, gtruth);

#if _DLEVEL > 2
	gtruth_print_vehicles(vt);
#endif

	return vt;
}

static void print_vehicle_info(const struct vehicle_mark *v)
{
	printf("   iframe: %05d\n", v->iframe);
	printf("     lane: %d\n", v->lane);
	if (v->lane) {
		printf("    plate: True\n");
	} else {
		printf("    plate: False\n");
	}
	printf("        x: %d\n", v->x);
	printf("        y: %d\n", v->y);
	printf("        w: %d\n", v->w);
	printf("        h: %d\n", v->h);
	if (v->radar) {
		printf("    radar: True\n");
		printf("    start: %05d\n", v->frame_start);
		printf("      end: %05d\n", v->frame_end);
		printf("    speed: %.2lf\n", v->speed);
	} else {
		printf("    radar: False\n");
	}
}

void gtruth_print_vehicles(const struct vehicle_table *vt)
{
	int i;

	printf("Print information about vehicle plates and groundtruth:\n\n");
	for (i = 0; i < vt->size; i++) {
		if (!vt->vehicles[i]) {
			continue;
		}
		print_vehicle_info(vt->vehicles[i]);
		printf("\n");
	}
	printf("\n------------------------------------------------------\n");
}


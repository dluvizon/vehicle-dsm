/*
 * GStreamer
 * Copyright (C) 2005 Thomas Vander Stichele <thomas@apestaart.org>
 * Copyright (C) 2005 Ronald S. Bultje <rbultje@ronald.bitfreak.net>
 * Copyright (C) YEAR AUTHOR_NAME AUTHOR_EMAIL
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * Alternatively, the contents of this file may be used under the
 * GNU Lesser General Public License Version 2.1 (the "LGPL"), in
 * which case the following provisions apply instead of the ones
 * mentioned above:
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/**
 * SECTION:element
 *
 * FIXME:Describe plugin here.
 *
 * <refsect2>
 * <title>Example launch line</title>
 * |[
 * gst-launch -v -m fakesrc ! vehicledsm ! fakesink silent=TRUE
 * ]|
 * </refsect2>
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "gst-vehicle-dsm.h"

#ifdef  __cplusplus
extern "C" void *vs_preload(struct dsm_common_data *cdata, int width,
		int height);
extern "C" void vs_chain_callback(struct dsm_common_data *cd, void *_data);
extern "C" void vs_end_of_stream_handler(void *_data);
#else
void *vs_preload(struct dsm_common_data *cdata, int width, int height);
void vs_chain_callback(struct dsm_common_data *cd, void *_data);
void vs_end_of_stream_handler(void *_data);
#endif

GST_DEBUG_CATEGORY_STATIC(gst_vehicle_dsm_debug);
#define GST_CAT_DEFAULT gst_vehicle_dsm_debug

#define GST_PROP_GT_PATH_DEFAULT ""
#define GST_PROP_VEHICLES_PATH_DEFAULT ""
#define GST_PROP_OUT_PATH_DEFAULT "."

/* Filter signals and args */
enum {
	LAST_SIGNAL
};

enum {
	PROP_0,
	PROP_SILENT,
	PROP_GT_PATH,
	PROP_VEHICLES_PATH,
	PROP_OUT_PATH,
};

/* the capabilities of the inputs and outputs.
 * describe the real formats here.
 */
static GstStaticPadTemplate sink_factory = GST_STATIC_PAD_TEMPLATE(
		"sink",
		GST_PAD_SINK,
		GST_PAD_ALWAYS,
		GST_STATIC_CAPS ("ANY"));

static GstStaticPadTemplate src_factory = GST_STATIC_PAD_TEMPLATE(
		"src",
		GST_PAD_SRC,
		GST_PAD_ALWAYS,
		GST_STATIC_CAPS ("ANY"));

#define gst_vehicle_dsm_parent_class parent_class
G_DEFINE_TYPE(GstPluginVehicleDsm, gst_vehicle_dsm, GST_TYPE_ELEMENT);

static void gst_vehicle_dsm_set_property(
		GObject * object, guint prop_id,
		const GValue * value, GParamSpec * pspec);
static void gst_vehicle_dsm_get_property(
		GObject * object, guint prop_id,
		GValue * value, GParamSpec * pspec);

static gboolean gst_vehicle_dsm_sink_event(
		GstPad * pad, GstObject * parent, GstEvent * event);
static GstFlowReturn gst_vehicle_dsm_chain(
		GstPad * pad, GstObject * parent, GstBuffer * buf);

/* GObject vmethod implementations */

/* initialize the vehicledsm's class */
static void gst_vehicle_dsm_class_init(GstPluginVehicleDsmClass * klass)
{
	GObjectClass *gobject_class;
	GstElementClass *gstelement_class;

	gobject_class = (GObjectClass *) klass;
	gstelement_class = (GstElementClass *) klass;

	gobject_class->set_property = gst_vehicle_dsm_set_property;
	gobject_class->get_property = gst_vehicle_dsm_get_property;

	g_object_class_install_property (gobject_class, PROP_SILENT,
			g_param_spec_boolean (
				"silent", "Silent", "Produce verbose output",
				FALSE, G_PARAM_READWRITE));
	g_object_class_install_property (gobject_class, PROP_GT_PATH,
			g_param_spec_string (
				"groundtruth-path", "Ground Truth Path",
				"XML ground truth file",
				GST_PROP_GT_PATH_DEFAULT,
				G_PARAM_WRITABLE));
	g_object_class_install_property (gobject_class, PROP_VEHICLES_PATH,
			g_param_spec_string (
				"vehicles.xml", "Vehicles XML file",
				"XML file with plates and ground truth",
				GST_PROP_VEHICLES_PATH_DEFAULT,
				G_PARAM_WRITABLE));
	g_object_class_install_property (gobject_class, PROP_OUT_PATH,
			g_param_spec_string (
				"output-path", "Output Path",
				"Path to the output folder",
				GST_PROP_OUT_PATH_DEFAULT,
				G_PARAM_WRITABLE));

	gst_element_class_set_details_simple(gstelement_class,
			"Vehicle Detection and Speed Measurement",
			"Computer Vision",
			"This software can detect vehicles and estimate "
			"their speed",
			"Diogo Luvizon <diogo@luvizon.com>\n"
"                           Rodrigo Minetto <rodrigo.minetto@gmail.com>\n"
"                           Bogdan Nassu <btnassu@yahoo.com.br>");

	gst_element_class_add_pad_template (gstelement_class,
			gst_static_pad_template_get (&src_factory));
	gst_element_class_add_pad_template (gstelement_class,
			gst_static_pad_template_get (&sink_factory));
}

/* initialize the new element
 * instantiate pads and add them to element
 * set pad calback functions
 * initialize instance structure
 */
static void gst_vehicle_dsm_init(GstPluginVehicleDsm *filter)
{
	filter->sinkpad = gst_pad_new_from_static_template(&sink_factory,
			"sink");
	gst_pad_set_event_function (filter->sinkpad,
			GST_DEBUG_FUNCPTR(gst_vehicle_dsm_sink_event));
	gst_pad_set_chain_function (filter->sinkpad,
			GST_DEBUG_FUNCPTR(gst_vehicle_dsm_chain));
	GST_PAD_SET_PROXY_CAPS (filter->sinkpad);
	gst_element_add_pad (GST_ELEMENT (filter), filter->sinkpad);

	filter->srcpad = gst_pad_new_from_static_template (&src_factory,
			"src");
	GST_PAD_SET_PROXY_CAPS (filter->srcpad);
	gst_element_add_pad (GST_ELEMENT (filter), filter->srcpad);

	filter->silent = FALSE;
	filter->cdata.groundtruh_file = GST_PROP_GT_PATH_DEFAULT;
	filter->cdata.vehicles_xml = GST_PROP_VEHICLES_PATH_DEFAULT;
	filter->cdata.output_path = GST_PROP_OUT_PATH_DEFAULT;
	filter->cdata.iframe = 0;
	filter->cdata.icurr = 0;
	filter->cdata.nhold = 0;
}

static void gst_vehicle_dsm_set_property(GObject * object,
		guint prop_id, const GValue * value, GParamSpec * pspec)
{
	size_t len;
	gchar *ps;
	GstPluginVehicleDsm *filter = GST_VEHICLE_DSM (object);

	switch (prop_id) {
	case PROP_SILENT:
		filter->silent = g_value_get_boolean (value);
		break;
	case PROP_GT_PATH:
		ps = g_value_dup_string(value);
		len = strlen(ps);
		if (ps[len - 1] == '/') {
			/* remove trailing slash if exist */
			ps = '\0';
		}
		filter->cdata.groundtruh_file = ps;
		break;
	case PROP_VEHICLES_PATH:
		ps = g_value_dup_string(value);
		len = strlen(ps);
		if (ps[len - 1] == '/') {
			/* remove trailing slash if exist */
			ps = '\0';
		}
		filter->cdata.vehicles_xml = ps;
		break;
	case PROP_OUT_PATH:
		ps = g_value_dup_string(value);
		len = strlen(ps);
		if (ps[len - 1] == '/') {
			/* remove trailing slash if exist */
			ps[len - 1] = '\0';
		}
		filter->cdata.output_path = ps;
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	break;
	}
}

static void gst_vehicle_dsm_get_property(GObject * object,
		guint prop_id, GValue * value, GParamSpec * pspec)
{
	GstPluginVehicleDsm *filter = GST_VEHICLE_DSM (object);

	switch (prop_id) {
	case PROP_SILENT:
		g_value_set_boolean (value, filter->silent);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

/* GstElement vmethod implementations */

/* this function handles sink events */
static gboolean gst_vehicle_dsm_sink_event(GstPad *pad,
		GstObject *parent, GstEvent *event)
{
	gboolean ret;
	GstPluginVehicleDsm *filter;

	filter = GST_VEHICLE_DSM(parent);

	switch (GST_EVENT_TYPE(event)) {
	case GST_EVENT_CAPS: {
		GstCaps * caps;

		gst_event_parse_caps (event, &caps);
		/* do something with the caps */

		/* and forward */
		ret = gst_pad_event_default(pad, parent, event);
		break;
	} case GST_EVENT_EOS: {
		g_print("Received EOS event!\n");
		vs_end_of_stream_handler(filter->_priv_data);
		ret = TRUE;
		break;
	} default:
		ret = gst_pad_event_default(pad, parent, event);
		break;
	}
	return ret;
}

/* chain function
 * this function does the actual processing
 */
static GstFlowReturn gst_vehicle_dsm_chain(GstPad *pad,
		GstObject *parent, GstBuffer *buf)
{
	GstPluginVehicleDsm *filter;
	const GstStructure *str;
	struct dsm_common_data *cd;
	GstFlowReturn ret;
	GstCaps *caps;
	GstMapInfo *info;
	gint i;

	ret = GST_FLOW_CUSTOM_SUCCESS;

	/* retrive data from GstObject */
	filter = GST_VEHICLE_DSM(parent);
	cd = &filter->cdata;

	i = cd->icurr;
	filter->buf[i] = buf;
	info = &filter->info[i];
	cd->nhold++;
	gst_buffer_map(buf, info, GST_MAP_READ);
	cd->vdata[i] = info->data;

	if (0 == cd->iframe) {
		caps = gst_pad_get_current_caps(pad);
		str = gst_caps_get_structure(caps, 0);
		if (!gst_structure_get_int(str, "width", &filter->width) ||
				!gst_structure_get_int(str, "height",
					&filter->height)) {
			g_error("No width/height available\n");
		}
		filter->_priv_data = (gpointer) vs_preload(&filter->cdata,
				filter->width, filter->height);
		if (!filter->_priv_data) {
			exit(1);
		}
	}

	if (cd->nhold >= NUM_OF_STILL_BUF) {
		vs_chain_callback(cd, filter->_priv_data);
	}
	cd->iframe++;
	cd->icurr++;
	if (cd->icurr >= NUM_OF_STILL_BUF) {
		cd->icurr %= NUM_OF_STILL_BUF;
	}
	if (cd->nhold >= NUM_OF_STILL_BUF) {
		/* if we are filled, start calling the next chain function */
		i = cd->icurr;
		info = &filter->info[i];

		gst_buffer_unmap(filter->buf[i], info);
		cd->nhold--;

		ret = gst_pad_push(filter->srcpad, filter->buf[i]);
	}


	return ret;
}

/* entry point to initialize the plug-in
 * initialize the plug-in itself
 * register the element factories and other features
 */
static gboolean vehicle_dsm_init (GstPlugin * plugin)
{
	/* debug category for fltering log messages
	 *
	 * exchange the string 'VehicleDsm plugin' with your description
	 */
	GST_DEBUG_CATEGORY_INIT (gst_vehicle_dsm_debug,
			"vehicledsm", 0, "Vehicle DSM plugin");

	return gst_element_register(plugin, "vehicledsm",
			GST_RANK_NONE, GST_TYPE_VEHICLE_DSM);
}

/* PACKAGE: this is usually set by autotools depending on some _INIT macro
 * in configure.ac and then written into and defined in config.h, but we can
 * just set it ourselves here in case someone doesn't use autotools to
 * compile this code. GST_DEFINE needs PACKAGE to be defined.
 */
#ifndef PACKAGE
#define PACKAGE "vehicledsm"
#endif

/* gstreamer looks for this structure to register plugins
 *
 * exchange the string 'VehicleDsm plugin' with your plugin description
 */
GST_PLUGIN_DEFINE(
		GST_VERSION_MAJOR,
		GST_VERSION_MINOR,
		vehicledsm,
		"Vehicle DSM plugin",
		vehicle_dsm_init,
		VERSION,
		"LGPL",
		"GStreamer",
		"http://gstreamer.net/");


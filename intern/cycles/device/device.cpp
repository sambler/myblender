/*
 * Copyright 2011-2013 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>
#include <string.h>

#include "device.h"
#include "device_intern.h"

#include "util_debug.h"
#include "util_foreach.h"
#include "util_half.h"
#include "util_math.h"
#include "util_opengl.h"
#include "util_time.h"
#include "util_types.h"
#include "util_vector.h"

CCL_NAMESPACE_BEGIN

/* Device Requested Features */

std::ostream& operator <<(std::ostream &os,
                          const DeviceRequestedFeatures& requested_features)
{
	os << "Experimental features: "
	   << (requested_features.experimental ? "On" : "Off") << std::endl;
	os << "Max closure count: " << requested_features.max_closure << std::endl;
	os << "Max nodes group: " << requested_features.max_nodes_group << std::endl;
	/* TODO(sergey): Decode bitflag into list of names. */
	os << "Nodes features: " << requested_features.nodes_features << std::endl;
	/* TODO(sergey): Make it utility function to convert bool to string. */
	os << "Use hair: "
	   << (requested_features.use_hair ? "True" : "False")  << std::endl;
	os << "Use object motion: "
	   << (requested_features.use_object_motion ? "True" : "False")  << std::endl;
	os << "Use camera motion: "
	   << (requested_features.use_camera_motion ? "True" : "False")  << std::endl;
	os << "Use Baking: "
	   << (requested_features.use_baking ? "True" : "False")  << std::endl;
	return os;
}

/* Device */

Device::~Device()
{
	if(!background && vertex_buffer != 0) {
		glDeleteBuffers(1, &vertex_buffer);
	}
}

void Device::pixels_alloc(device_memory& mem)
{
	mem_alloc(mem, MEM_READ_WRITE);
}

void Device::pixels_copy_from(device_memory& mem, int y, int w, int h)
{
	if(mem.data_type == TYPE_HALF)
		mem_copy_from(mem, y, w, h, sizeof(half4));
	else
		mem_copy_from(mem, y, w, h, sizeof(uchar4));
}

void Device::pixels_free(device_memory& mem)
{
	mem_free(mem);
}

void Device::draw_pixels(device_memory& rgba, int y, int w, int h, int dx, int dy, int width, int height, bool transparent,
	const DeviceDrawParams &draw_params)
{
	pixels_copy_from(rgba, y, w, h);

	if(transparent) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	}

	glColor3f(1.0f, 1.0f, 1.0f);

	if(rgba.data_type == TYPE_HALF) {
		/* for multi devices, this assumes the inefficient method that we allocate
		 * all pixels on the device even though we only render to a subset */
		GLhalf *data_pointer = (GLhalf*)rgba.data_pointer;
		float vbuffer[16], *basep;
		float *vp = NULL;

		data_pointer += 4*y*w;

		/* draw half float texture, GLSL shader for display transform assumed to be bound */
		GLuint texid;
		glGenTextures(1, &texid);
		glBindTexture(GL_TEXTURE_2D, texid);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, w, h, 0, GL_RGBA, GL_HALF_FLOAT, data_pointer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		glEnable(GL_TEXTURE_2D);

		if(draw_params.bind_display_space_shader_cb) {
			draw_params.bind_display_space_shader_cb();
		}

		if(GLEW_VERSION_1_5) {
			if(!vertex_buffer)
				glGenBuffers(1, &vertex_buffer);

			glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
			/* invalidate old contents - avoids stalling if buffer is still waiting in queue to be rendered */
			glBufferData(GL_ARRAY_BUFFER, 16 * sizeof(float), NULL, GL_STREAM_DRAW);

			vp = (float *)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

			basep = NULL;
		}
		else {
			basep = vbuffer;
			vp = vbuffer;
		}

		if(vp) {
			/* texture coordinate - vertex pair */
			vp[0] = 0.0f;
			vp[1] = 0.0f;
			vp[2] = dx;
			vp[3] = dy;

			vp[4] = 1.0f;
			vp[5] = 0.0f;
			vp[6] = (float)width + dx;
			vp[7] = dy;

			vp[8] = 1.0f;
			vp[9] = 1.0f;
			vp[10] = (float)width + dx;
			vp[11] = (float)height + dy;

			vp[12] = 0.0f;
			vp[13] = 1.0f;
			vp[14] = dx;
			vp[15] = (float)height + dy;

			if(vertex_buffer)
				glUnmapBuffer(GL_ARRAY_BUFFER);
		}

		glTexCoordPointer(2, GL_FLOAT, 4 * sizeof(float), basep);
		glVertexPointer(2, GL_FLOAT, 4 * sizeof(float), ((char *)basep) + 2 * sizeof(float));

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);

		glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);

		if(vertex_buffer) {
			glBindBuffer(GL_ARRAY_BUFFER, 0);
		}

		if(draw_params.unbind_display_space_shader_cb) {
			draw_params.unbind_display_space_shader_cb();
		}

		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_TEXTURE_2D);
		glDeleteTextures(1, &texid);
	}
	else {
		/* fallback for old graphics cards that don't support GLSL, half float,
		 * and non-power-of-two textures */
		glPixelZoom((float)width/(float)w, (float)height/(float)h);
		glRasterPos2f(dx, dy);

		uint8_t *pixels = (uint8_t*)rgba.data_pointer;

		pixels += 4*y*w;

		glDrawPixels(w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

		glRasterPos2f(0.0f, 0.0f);
		glPixelZoom(1.0f, 1.0f);
	}

	if(transparent)
		glDisable(GL_BLEND);
}

Device *Device::create(DeviceInfo& info, Stats &stats, bool background)
{
	Device *device;

	switch(info.type) {
		case DEVICE_CPU:
			device = device_cpu_create(info, stats, background);
			break;
#ifdef WITH_CUDA
		case DEVICE_CUDA:
			if(device_cuda_init())
				device = device_cuda_create(info, stats, background);
			else
				device = NULL;
			break;
#endif
#ifdef WITH_MULTI
		case DEVICE_MULTI:
			device = device_multi_create(info, stats, background);
			break;
#endif
#ifdef WITH_NETWORK
		case DEVICE_NETWORK:
			device = device_network_create(info, stats, "127.0.0.1");
			break;
#endif
#ifdef WITH_OPENCL
		case DEVICE_OPENCL:
			if(device_opencl_init())
				device = device_opencl_create(info, stats, background);
			else
				device = NULL;
			break;
#endif
		default:
			return NULL;
	}

	return device;
}

DeviceType Device::type_from_string(const char *name)
{
	if(strcmp(name, "cpu") == 0)
		return DEVICE_CPU;
	else if(strcmp(name, "cuda") == 0)
		return DEVICE_CUDA;
	else if(strcmp(name, "opencl") == 0)
		return DEVICE_OPENCL;
	else if(strcmp(name, "network") == 0)
		return DEVICE_NETWORK;
	else if(strcmp(name, "multi") == 0)
		return DEVICE_MULTI;
	
	return DEVICE_NONE;
}

string Device::string_from_type(DeviceType type)
{
	if(type == DEVICE_CPU)
		return "cpu";
	else if(type == DEVICE_CUDA)
		return "cuda";
	else if(type == DEVICE_OPENCL)
		return "opencl";
	else if(type == DEVICE_NETWORK)
		return "network";
	else if(type == DEVICE_MULTI)
		return "multi";
	
	return "";
}

vector<DeviceType>& Device::available_types()
{
	static vector<DeviceType> types;
	static bool types_init = false;

	if(!types_init) {
		types.push_back(DEVICE_CPU);

#ifdef WITH_CUDA
		if(device_cuda_init())
			types.push_back(DEVICE_CUDA);
#endif

#ifdef WITH_OPENCL
		if(device_opencl_init())
			types.push_back(DEVICE_OPENCL);
#endif

#ifdef WITH_NETWORK
		types.push_back(DEVICE_NETWORK);
#endif
#ifdef WITH_MULTI
		types.push_back(DEVICE_MULTI);
#endif

		types_init = true;
	}

	return types;
}

vector<DeviceInfo>& Device::available_devices()
{
	static vector<DeviceInfo> devices;
	static bool devices_init = false;

	if(!devices_init) {
#ifdef WITH_CUDA
		if(device_cuda_init())
			device_cuda_info(devices);
#endif

#ifdef WITH_OPENCL
		if(device_opencl_init())
			device_opencl_info(devices);
#endif

#ifdef WITH_MULTI
		device_multi_info(devices);
#endif

		device_cpu_info(devices);

#ifdef WITH_NETWORK
		device_network_info(devices);
#endif

		devices_init = true;
	}

	return devices;
}

string Device::device_capabilities()
{
	string capabilities = "CPU device capabilities: ";
	capabilities += device_cpu_capabilities() + "\n";
#ifdef WITH_CUDA
	if(device_cuda_init()) {
		capabilities += "\nCUDA device capabilities:\n";
		capabilities += device_cuda_capabilities();
	}
#endif

#ifdef WITH_OPENCL
	if(device_opencl_init()) {
		capabilities += "\nOpenCL device capabilities:\n";
		capabilities += device_opencl_capabilities();
	}
#endif

	return capabilities;
}

CCL_NAMESPACE_END

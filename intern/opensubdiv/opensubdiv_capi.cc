/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2013 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Sergey Sharybin.
 *                 Brecht van Lommel
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include "opensubdiv_capi.h"

#ifdef _MSC_VER
#  include "iso646.h"
#endif

#include <stdlib.h>
#include <GL/glew.h>

#include <opensubdiv/osd/glMesh.h>

/* CPU Backend */
#include <opensubdiv/osd/cpuGLVertexBuffer.h>
#include <opensubdiv/osd/cpuEvaluator.h>

#ifdef OPENSUBDIV_HAS_OPENMP
#  include <opensubdiv/osd/ompEvaluator.h>
#endif  /* OPENSUBDIV_HAS_OPENMP */

#ifdef OPENSUBDIV_HAS_OPENCL
#  include <opensubdiv/osd/clGLVertexBuffer.h>
#  include <opensubdiv/osd/clEvaluator.h>
#  include "opensubdiv_device_context_opencl.h"
#endif  /* OPENSUBDIV_HAS_OPENCL */

#ifdef OPENSUBDIV_HAS_CUDA
#  include <opensubdiv/osd/cudaGLVertexBuffer.h>
#  include <opensubdiv/osd/cudaEvaluator.h>
#  include "opensubdiv_device_context_cuda.h"
#endif  /* OPENSUBDIV_HAS_CUDA */

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
#  include <opensubdiv/osd/glXFBEvaluator.h>
#  include <opensubdiv/osd/glVertexBuffer.h>
#endif  /* OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK */

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
#  include <opensubdiv/osd/glComputeEvaluator.h>
#  include <opensubdiv/osd/glVertexBuffer.h>
#endif  /* OPENSUBDIV_HAS_GLSL_COMPUTE */

#include <opensubdiv/osd/glPatchTable.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/primvarRefiner.h>

#include "opensubdiv_intern.h"
#include "opensubdiv_topology_refiner.h"

#include "MEM_guardedalloc.h"

/* **************** Types declaration **************** */

using OpenSubdiv::Osd::GLMeshInterface;
using OpenSubdiv::Osd::Mesh;
using OpenSubdiv::Osd::MeshBitset;
using OpenSubdiv::Far::StencilTable;
using OpenSubdiv::Osd::GLPatchTable;

using OpenSubdiv::Osd::Mesh;

/* CPU backend */
using OpenSubdiv::Osd::CpuGLVertexBuffer;
using OpenSubdiv::Osd::CpuEvaluator;
typedef Mesh<CpuGLVertexBuffer,
             StencilTable,
             CpuEvaluator,
             GLPatchTable> OsdCpuMesh;

#ifdef OPENSUBDIV_HAS_OPENMP
using OpenSubdiv::Osd::OmpEvaluator;
typedef Mesh<CpuGLVertexBuffer,
             StencilTable,
             OmpEvaluator,
             GLPatchTable> OsdOmpMesh;
#endif  /* OPENSUBDIV_HAS_OPENMP */

#ifdef OPENSUBDIV_HAS_OPENCL
using OpenSubdiv::Osd::CLEvaluator;
using OpenSubdiv::Osd::CLGLVertexBuffer;
using OpenSubdiv::Osd::CLStencilTable;
/* TODO(sergey): Use CLDeviceCOntext similar to OSD examples? */
typedef Mesh<CLGLVertexBuffer,
             CLStencilTable,
             CLEvaluator,
             GLPatchTable,
             CLDeviceContext> OsdCLMesh;
static CLDeviceContext g_clDeviceContext;
#endif  /* OPENSUBDIV_HAS_OPENCL */

#ifdef OPENSUBDIV_HAS_CUDA
using OpenSubdiv::Osd::CudaEvaluator;
using OpenSubdiv::Osd::CudaGLVertexBuffer;
using OpenSubdiv::Osd::CudaStencilTable;
typedef Mesh<CudaGLVertexBuffer,
             CudaStencilTable,
             CudaEvaluator,
             GLPatchTable> OsdCudaMesh;
static CudaDeviceContext g_cudaDeviceContext;
#endif  /* OPENSUBDIV_HAS_CUDA */

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
using OpenSubdiv::Osd::GLXFBEvaluator;
using OpenSubdiv::Osd::GLStencilTableTBO;
using OpenSubdiv::Osd::GLVertexBuffer;
typedef Mesh<GLVertexBuffer,
             GLStencilTableTBO,
             GLXFBEvaluator,
             GLPatchTable> OsdGLSLTransformFeedbackMesh;
#endif  /* OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK */

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
using OpenSubdiv::Osd::GLComputeEvaluator;
using OpenSubdiv::Osd::GLStencilTableSSBO;
using OpenSubdiv::Osd::GLVertexBuffer;
typedef Mesh<GLVertexBuffer,
             GLStencilTableSSBO,
             GLComputeEvaluator,
             GLPatchTable> OsdGLSLComputeMesh;
#endif

namespace {

struct FVarVertex {
	float u, v;
	void Clear() {
		u = v = 0.0f;
	}
	void AddWithWeight(FVarVertex const & src, float weight) {
		u += weight * src.u;
		v += weight * src.v;
	}
};

static void interpolate_fvar_data(OpenSubdiv::Far::TopologyRefiner& refiner,
                                  const std::vector<float> uvs,
                                  std::vector<float> &fvar_data) {
	/* TODO(sergey): Support all FVar channels. */
	const int channel = 0;
	/* TODO(sergey): Make it somehow more generic way. */
	const int fvar_width = 2;

	int max_level = refiner.GetMaxLevel(),
	    num_values_max = refiner.GetLevel(max_level).GetNumFVarValues(channel),
	    num_values_total = refiner.GetNumFVarValuesTotal(channel);
	if (num_values_total <= 0) {
		return;
	}
	OpenSubdiv::Far::PrimvarRefiner primvarRefiner(refiner);
	if (refiner.IsUniform()) {
		/* For uniform we only keep the highest level of refinement. */
		fvar_data.resize(num_values_max * fvar_width);
		std::vector<FVarVertex> buffer(num_values_total - num_values_max);
		FVarVertex *src = &buffer[0];
		memcpy(src, &uvs[0], uvs.size() * sizeof(float));
		/* Defer the last level to treat separately with its alternate
		 * destination.
		 */
		for (int level = 1; level < max_level; ++level) {
			FVarVertex *dst = src + refiner.GetLevel(level-1).GetNumFVarValues(channel);
			primvarRefiner.InterpolateFaceVarying(level, src, dst, channel);
			src = dst;
		}
		FVarVertex *dst = reinterpret_cast<FVarVertex *>(&fvar_data[0]);
		primvarRefiner.InterpolateFaceVarying(max_level, src, dst, channel);
	} else {
		/* For adaptive we keep all levels. */
		fvar_data.resize(num_values_total * fvar_width);
		FVarVertex *src = reinterpret_cast<FVarVertex *>(&fvar_data[0]);
		memcpy(src, &uvs[0], uvs.size() * sizeof(float));
		for (int level = 1; level <= max_level; ++level) {
			FVarVertex *dst = src + refiner.GetLevel(level-1).GetNumFVarValues(channel);
			primvarRefiner.InterpolateFaceVarying(level, src, dst, channel);
			src = dst;
        }
    }
}

}  // namespace

struct OpenSubdiv_GLMesh *openSubdiv_createOsdGLMeshFromTopologyRefiner(
        OpenSubdiv_TopologyRefinerDescr *topology_refiner,
        int evaluator_type,
        int level)
{
	using OpenSubdiv::Far::TopologyRefiner;

	MeshBitset bits;
	/* TODO(sergey): Adaptive subdivisions are not currently
	 * possible because of the lack of tessellation shader.
	 */
	bits.set(OpenSubdiv::Osd::MeshAdaptive, 0);
	bits.set(OpenSubdiv::Osd::MeshUseSingleCreasePatch, 0);
	bits.set(OpenSubdiv::Osd::MeshInterleaveVarying, 1);
	bits.set(OpenSubdiv::Osd::MeshFVarData, 1);
	bits.set(OpenSubdiv::Osd::MeshEndCapBSplineBasis, 1);

	const int num_vertex_elements = 3;
	const int num_varying_elements = 3;

	GLMeshInterface *mesh = NULL;
	TopologyRefiner *refiner = topology_refiner->osd_refiner;

	switch(evaluator_type) {
#define CHECK_EVALUATOR_TYPE(type, class) \
		case OPENSUBDIV_EVALUATOR_ ## type: \
			mesh = new class(refiner, \
			                 num_vertex_elements, \
			                 num_varying_elements, \
			                 level, \
			                 bits); \
			break;

		CHECK_EVALUATOR_TYPE(CPU, OsdCpuMesh)

#ifdef OPENSUBDIV_HAS_OPENMP
		CHECK_EVALUATOR_TYPE(OPENMP, OsdOmpMesh)
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
		CHECK_EVALUATOR_TYPE(OPENCL, OsdCLMesh)
#endif

#ifdef OPENSUBDIV_HAS_CUDA
		CHECK_EVALUATOR_TYPE(CUDA, OsdCudaMesh)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
		CHECK_EVALUATOR_TYPE(GLSL_TRANSFORM_FEEDBACK,
		                     OsdGLSLTransformFeedbackMesh)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
		CHECK_EVALUATOR_TYPE(GLSL_COMPUTE, OsdGLSLComputeMesh)
#endif

#undef CHECK_EVALUATOR_TYPE
	}

	if (mesh == NULL) {
		return NULL;
	}

	OpenSubdiv_GLMesh *gl_mesh =
		(OpenSubdiv_GLMesh *) OBJECT_GUARDED_NEW(OpenSubdiv_GLMesh);
	gl_mesh->evaluator_type = evaluator_type;
	gl_mesh->descriptor = (OpenSubdiv_GLMeshDescr *) mesh;
	gl_mesh->topology_refiner = topology_refiner;

	if (refiner->GetNumFVarChannels() > 0) {
		std::vector<float> fvar_data;
		interpolate_fvar_data(*refiner, topology_refiner->uvs, fvar_data);
		openSubdiv_osdGLAllocFVar(gl_mesh, &fvar_data[0]);
	}
	else {
		gl_mesh->fvar_data = NULL;
	}

	return gl_mesh;
}

void openSubdiv_deleteOsdGLMesh(struct OpenSubdiv_GLMesh *gl_mesh)
{
	openSubdiv_osdGLDestroyFVar(gl_mesh);
	switch (gl_mesh->evaluator_type) {
#define CHECK_EVALUATOR_TYPE(type, class) \
		case OPENSUBDIV_EVALUATOR_ ## type: \
			delete (class *) gl_mesh->descriptor; \
			break;

		CHECK_EVALUATOR_TYPE(CPU, OsdCpuMesh)

#ifdef OPENSUBDIV_HAS_OPENMP
		CHECK_EVALUATOR_TYPE(OPENMP, OsdOmpMesh)
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
		CHECK_EVALUATOR_TYPE(OPENCL, OsdCLMesh)
#endif

#ifdef OPENSUBDIV_HAS_CUDA
		CHECK_EVALUATOR_TYPE(CUDA, OsdCudaMesh)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
		CHECK_EVALUATOR_TYPE(GLSL_TRANSFORM_FEEDBACK,
		                     OsdGLSLTransformFeedbackMesh)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
		CHECK_EVALUATOR_TYPE(GLSL_COMPUTE, OsdGLSLComputeMesh)
#endif

#undef CHECK_EVALUATOR_TYPE
	}

	/* NOTE: OSD refiner was owned by gl_mesh, no need to free it here. */
	OBJECT_GUARDED_DELETE(gl_mesh->topology_refiner, OpenSubdiv_TopologyRefinerDescr);
	OBJECT_GUARDED_DELETE(gl_mesh, OpenSubdiv_GLMesh);
}

unsigned int openSubdiv_getOsdGLMeshPatchIndexBuffer(struct OpenSubdiv_GLMesh *gl_mesh)
{
	return ((GLMeshInterface *)gl_mesh->descriptor)->GetPatchTable()->GetPatchIndexBuffer();
}

unsigned int openSubdiv_getOsdGLMeshVertexBuffer(struct OpenSubdiv_GLMesh *gl_mesh)
{
	return ((GLMeshInterface *)gl_mesh->descriptor)->BindVertexBuffer();
}

void openSubdiv_osdGLMeshUpdateVertexBuffer(struct OpenSubdiv_GLMesh *gl_mesh,
                                            const float *vertex_data,
                                            int start_vertex,
                                            int num_verts)
{
	((GLMeshInterface *)gl_mesh->descriptor)->UpdateVertexBuffer(vertex_data,
	                                                             start_vertex,
	                                                             num_verts);
}

void openSubdiv_osdGLMeshRefine(struct OpenSubdiv_GLMesh *gl_mesh)
{
	((GLMeshInterface *)gl_mesh->descriptor)->Refine();
}

void openSubdiv_osdGLMeshSynchronize(struct OpenSubdiv_GLMesh *gl_mesh)
{
	((GLMeshInterface *)gl_mesh->descriptor)->Synchronize();
}

void openSubdiv_osdGLMeshBindVertexBuffer(OpenSubdiv_GLMesh *gl_mesh)
{
	((GLMeshInterface *)gl_mesh->descriptor)->BindVertexBuffer();
}

const struct OpenSubdiv_TopologyRefinerDescr *openSubdiv_getGLMeshTopologyRefiner(
        OpenSubdiv_GLMesh *gl_mesh)
{
	return gl_mesh->topology_refiner;
}

int openSubdiv_supportGPUDisplay(void)
{
	// TODO: simplify extension check once Blender adopts GL 3.2
	return openSubdiv_gpu_legacy_support() &&
	       (GLEW_VERSION_3_2 ||
	       (GLEW_VERSION_3_1 && GLEW_EXT_geometry_shader4) ||
	       (GLEW_VERSION_3_0 && GLEW_EXT_geometry_shader4 && GLEW_ARB_uniform_buffer_object && (GLEW_ARB_texture_buffer_object || GLEW_EXT_texture_buffer_object)));
	/* also ARB_explicit_attrib_location? */
}

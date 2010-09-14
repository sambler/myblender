/**
***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * Contributor(s): Miika Hämäläinen
 *
 * ***** END GPL LICENSE BLOCK *****
 */


#include "MEM_guardedalloc.h"

#include <math.h>
#include "stdio.h"

#include "BLI_blenlib.h"
#include "BLI_math.h"
#include "BLI_kdtree.h"

#include "BKE_blender.h"
#include "BKE_cdderivedmesh.h"
#include "BKE_context.h"
#include "BKE_customdata.h"
#include "BKE_colortools.h"
#include "BKE_DerivedMesh.h"
#include "BKE_global.h"
#include "BKE_modifier.h"
#include "BKE_particle.h"
#include "BKE_material.h"
#include "BKE_dynamicpaint.h"
#include "BKE_texture.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"
#include "DNA_dynamicpaint_types.h"

// to get temp path
#include "DNA_userdef_types.h"

// for bake operator
#include "ED_screen.h"
#include "WM_types.h"
#include "WM_api.h"

// for image output
#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"
#include "BKE_image.h"
#include "intern/IMB_filetype.h"
#ifdef WITH_OPENEXR
#include "intern/openexr/openexr_api.h"
#endif

// Platform independend time
#include "PIL_time.h"

// bvh tree
#include "BKE_bvhutils.h"

// to read object material color

#include "DNA_texture_types.h"
#include "../render/intern/include/render_types.h"
#include "DNA_material_types.h"
#include "RE_render_ext.h"


#define MPOUTPUT_JPEG 0
#define MPOUTPUT_PNG 1
#define MPOUTPUT_OPENEXR 2

#define MPOUTPUT_PAINT 0
#define MPOUTPUT_WET 1
#define MPOUTPUT_DISPLACE 2

struct Object;
struct Scene;
struct DerivedMesh;
struct DynamicPaintModifierData;

/*
*	Init predefined antialias jitter data
*/
float jitterDistances[5] = {0.0f,
							0.447213595f,
							0.447213595f,
							0.447213595f,
							0.5f};

// precalculated gaussian factors for 5x super sampling
float gaussianFactors[5] = {	0.996849f,
								0.596145f,
								0.596145f,
								0.596145f,
								0.524141f};
float gaussianTotal = 3.309425f;

/*
*	Neighbour pixel table x and y list
*/
int neighX[8] = {1,1,0,-1,-1,-1, 0, 1};
int neighY[8] = {0,1,1, 1, 0,-1,-1,-1};


/*
*	Update derived mesh data if baking
*/
void dynamicPaint_Modifier_update(DynamicPaintModifierData *pmd, DerivedMesh *dm) {

	if (!pmd->baking) return;

	if((pmd->type & MOD_DYNAMICPAINT_TYPE_CANVAS) && pmd->canvas) {

		if (pmd->canvas->dm) pmd->canvas->dm->release(pmd->canvas->dm);
		pmd->canvas->dm = CDDM_copy(dm);
	}
	else if((pmd->type & MOD_DYNAMICPAINT_TYPE_PAINT) && pmd->paint) {

		if (pmd->paint->dm) pmd->paint->dm->release(pmd->paint->dm);
		pmd->paint->dm = CDDM_copy(dm);
	}

}

/*
*	Free canvas data
*/
static void dynamicPaint_Modifier_freeCanvas(DynamicPaintModifierData *pmd)
{
	if(pmd->canvas)
	{
		if (pmd->canvas->dm)
			pmd->canvas->dm->release(pmd->canvas->dm);
		pmd->canvas->dm = NULL;
	}
}

/*
*	Free paint data
*/
static void dynamicPaint_Modifier_freePaint(DynamicPaintModifierData *pmd)
{
	if(pmd->paint)
	{
		if(pmd->paint->dm)
			pmd->paint->dm->release(pmd->paint->dm);
		pmd->paint->dm = NULL;


		if(pmd->paint->paint_ramp)
			 MEM_freeN(pmd->paint->paint_ramp);
		pmd->paint->paint_ramp = NULL;
	}
}

/*
*	Free whole modifier
*/
void dynamicPaint_Modifier_free (DynamicPaintModifierData *pmd)
{
	if(pmd)
	{
		dynamicPaint_Modifier_freeCanvas(pmd);
		dynamicPaint_Modifier_freePaint(pmd);
	}
}

/*
*	Initialize modifier data
*/
void dynamicPaint_Modifier_createType(struct DynamicPaintModifierData *pmd)
{
	if(pmd)
	{
		if(pmd->type & MOD_DYNAMICPAINT_TYPE_CANVAS)
		{
			if(pmd->canvas)
				dynamicPaint_Modifier_freeCanvas(pmd);

			pmd->canvas = MEM_callocN(sizeof(DynamicPaintCanvasSettings), "DynamicPaintCanvas");
			pmd->canvas->pmd = pmd;

			pmd->canvas->flags = MOD_DPAINT_ANTIALIAS | MOD_DPAINT_MULALPHA | MOD_DPAINT_DRY_LOG;
			pmd->canvas->output = MOD_DPAINT_OUT_PAINT;
			pmd->canvas->effect = 0;
			pmd->canvas->effect_ui = 1;

			pmd->canvas->diss_speed = 300;
			pmd->canvas->dry_speed = 300;
			pmd->canvas->dflat_speed = 300;

			pmd->canvas->disp_depth = 1.0f;
			pmd->canvas->disp_type = MOD_DPAINT_DISP_DISPLACE;
			pmd->canvas->disp_format = MOD_DPAINT_DISPFOR_PNG;

			pmd->canvas->resolution = 256;
			pmd->canvas->start_frame = 1;
			pmd->canvas->end_frame = 100;
			pmd->canvas->substeps = 0;

			pmd->canvas->spread_speed = 1.0f;
			pmd->canvas->drip_speed = 1.0f;
			pmd->canvas->shrink_speed = 1.0f;

			sprintf(pmd->canvas->paint_output_path, "%spaintmap", U.tempdir);
			sprintf(pmd->canvas->wet_output_path, "%swetmap", U.tempdir);
			sprintf(pmd->canvas->displace_output_path, "%sdispmap", U.tempdir);


			pmd->canvas->dm = NULL;

		}
		else if(pmd->type & MOD_DYNAMICPAINT_TYPE_PAINT)
		{
			if(pmd->paint)
				dynamicPaint_Modifier_freePaint(pmd);

			pmd->paint = MEM_callocN(sizeof(DynamicPaintPainterSettings), "DynamicPaint Paint");
			pmd->paint->pmd = pmd;

			pmd->paint->psys = NULL;

			pmd->paint->flags = MOD_DPAINT_DO_PAINT | MOD_DPAINT_DO_WETNESS | MOD_DPAINT_DO_DISPLACE;
			pmd->paint->collision = MOD_DPAINT_COL_VOLUME;
			
			pmd->paint->r = 1.0f;
			pmd->paint->g = 1.0f;
			pmd->paint->b = 1.0f;
			pmd->paint->alpha = 1.0f;
			pmd->paint->wetness = 1.0f;

			pmd->paint->paint_distance = 0.1f;
			pmd->paint->proximity_falloff = MOD_DPAINT_PRFALL_SHARP;

			pmd->paint->displace_distance = 0.5f;
			pmd->paint->prox_displace_strength = 0.5f;

			pmd->paint->particle_radius = 0.2;
			pmd->paint->particle_smooth = 0.05;

			pmd->paint->dm = NULL;

			/*
			*	Paint proximity falloff colorramp.
			*/
			{
				CBData *ramp;

				pmd->paint->paint_ramp = add_colorband(0);
				ramp = pmd->paint->paint_ramp->data;
				// Add default smooth-falloff ramp
				ramp[0].r = ramp[0].g = ramp[0].b = ramp[0].a = 1.0f;
				ramp[0].pos = 0.0f;
				ramp[1].r = ramp[1].g = ramp[1].b = ramp[1].pos = 1.0f;
				ramp[1].a = 0.0f;
				pmd->paint->paint_ramp->tot = 2;
			}


			/*
			*	Displace proximity edge-falloff curve.
			*/
			{
				pmd->paint->disp_curve = curvemapping_add(1, 0.0f, 1.0f, 1.0f, 0.2f);

				curvemapping_initialize(pmd->paint->disp_curve);
			}
		}
	}
}

void dynamicPaint_Modifier_copy(struct DynamicPaintModifierData *pmd, struct DynamicPaintModifierData *tpmd)
{
	// Init modifier
	tpmd->type = pmd->type;
	dynamicPaint_Modifier_createType(tpmd);

	// Copy data
	if (tpmd->canvas) {
		pmd->canvas->pmd = tpmd;

		tpmd->canvas->flags = pmd->canvas->flags;
		tpmd->canvas->output = pmd->canvas->output;
		tpmd->canvas->disp_type = pmd->canvas->disp_type;
		tpmd->canvas->disp_format = pmd->canvas->disp_format;
		tpmd->canvas->effect = pmd->canvas->effect;
		tpmd->canvas->effect_ui = 1;

		tpmd->canvas->resolution = pmd->canvas->resolution;
		tpmd->canvas->start_frame = pmd->canvas->start_frame;
		tpmd->canvas->end_frame = pmd->canvas->end_frame;
		tpmd->canvas->substeps = pmd->canvas->substeps;

		tpmd->canvas->dry_speed = pmd->canvas->dry_speed;
		tpmd->canvas->diss_speed = pmd->canvas->diss_speed;
		tpmd->canvas->disp_depth = pmd->canvas->disp_depth;
		tpmd->canvas->dflat_speed = pmd->canvas->dflat_speed;

		strncpy(tpmd->canvas->paint_output_path, pmd->canvas->paint_output_path, 240);
		strncpy(tpmd->canvas->wet_output_path, pmd->canvas->wet_output_path, 240);
		strncpy(tpmd->canvas->displace_output_path, pmd->canvas->displace_output_path, 240);

		tpmd->canvas->spread_speed = pmd->canvas->spread_speed;
		tpmd->canvas->drip_speed = pmd->canvas->drip_speed;
		tpmd->canvas->shrink_speed = pmd->canvas->shrink_speed;

	} else if (tpmd->paint) {
		pmd->canvas->pmd = tpmd;

		tpmd->paint->flags = pmd->paint->flags;
		tpmd->paint->collision = pmd->paint->collision;

		tpmd->paint->r = pmd->paint->r;
		tpmd->paint->g = pmd->paint->g;
		tpmd->paint->b = pmd->paint->b;
		tpmd->paint->alpha = pmd->paint->alpha;
		tpmd->paint->wetness = pmd->paint->wetness;

		tpmd->paint->particle_radius = pmd->paint->particle_radius;
		tpmd->paint->particle_smooth = pmd->paint->particle_smooth;
		tpmd->paint->paint_distance = pmd->paint->paint_distance;
		tpmd->paint->psys = pmd->paint->psys;
		tpmd->paint->displace_distance = pmd->paint->displace_distance;
		tpmd->paint->prox_displace_strength = pmd->paint->prox_displace_strength;

		tpmd->paint->paint_ramp = pmd->paint->paint_ramp;
		tpmd->paint->disp_curve = pmd->paint->disp_curve;

		tpmd->paint->proximity_falloff = pmd->paint->proximity_falloff;
	}
}



void dynamicPaint_Modifier_do(DynamicPaintModifierData *pmd, Scene *scene, Object *ob, DerivedMesh *dm)
{	
	// Update derived mesh data to modifier if baking
	dynamicPaint_Modifier_update(pmd, dm);
}

/* find closest point to p on line through l1,l2 and return lambda,
 * where (0 <= lambda <= 1) when cp is in the line segement l1,l2
 */
float closest_to_line_v2(float cp[2],float p[2], float l1[2], float l2[2])
{
	float h[3],u[3],lambda;
	sub_v2_v2v2(u, l2, l1);
	sub_v2_v2v2(h, p, l1);
	lambda =dot_v2v2(u,h)/dot_v2v2(u,u);
	cp[0] = l1[0] + u[0] * lambda;
	cp[1] = l1[1] + u[1] * lambda;
	return lambda;
}

/* point closest to v1 on line v2-v3 in 2D */
void closest_to_line_segment_v2(float *closest, float p[2], float l1[2], float l2[2])
{
	float lambda, cp[2];

	lambda= closest_to_line_v3(cp,p, l1, l2);

	if(lambda <= 0.0f)
		copy_v2_v2(closest, l1);
	else if(lambda >= 1.0f)
		copy_v2_v2(closest, l2);
	else
		copy_v2_v2(closest, cp);
}

int dynamicPaint_findNeighbourPixel(DynamicPaintCanvasSettings *canvas, int px, int py, int n_index) {
	// Only use face edges to detect neighbours. (Accurate enough and
	// faster/simplier than including possible tip directional faces)

	int numOfFaces = canvas->dm->getNumFaces(canvas->dm);
	MVert *mvert = NULL;
	MFace *mface = NULL;
	MTFace *tface = NULL;

	struct PaintSurface *surface;
	struct BB2d *faceBB = NULL;

	PaintSurfacePoint *cPoint;

	surface = canvas->surface;

	mvert = canvas->dm->getVertArray(canvas->dm);
	mface = canvas->dm->getFaceArray(canvas->dm);
	tface = DM_get_face_data_layer(canvas->dm, CD_MTFACE);



	/*
	*	Simple neighbour face finding algorithm:
	*	TODO: Implement something more accurate / optimized
	*/

	cPoint = (&surface->point[px+surface->w*py]);

	// Get closest edge to that subpixel on UV map
	{
		float pixel[2], dist, t_dist;
		int i, uindex[2], edge1_index, edge2_index, e1_index, e2_index, target_face;

		float closest_point[2], lambda, dir_vec[2];
		int target_uv1, target_uv2, final_pixel[2], final_index;

		float (*s_uv1),(*s_uv2), (*t_uv1), (*t_uv2);

		pixel[0] = ((float)(px + neighX[n_index]) + 0.5f) / (float)surface->w;
		pixel[1] = ((float)(py + neighY[n_index]) + 0.5f) / (float)surface->h;

		// Get uv indexes for current face part
		if (cPoint->quad) {
			uindex[0] = 0; uindex[1] = 2; uindex[2] = 3;
		}
		else {
			uindex[0] = 0; uindex[1] = 1; uindex[2] = 2;
		}

		/*
		*	Find closest edge to that pixel
		*/
		// Dist to first edge
		e1_index = cPoint->v1; e2_index = cPoint->v2; edge1_index = uindex[0]; edge2_index = uindex[1];
		dist = dist_to_line_segment_v2(pixel, tface[cPoint->index].uv[edge1_index], tface[cPoint->index].uv[edge2_index]);

		// Dist to second edge
		t_dist = dist_to_line_segment_v2(pixel, tface[cPoint->index].uv[uindex[1]], tface[cPoint->index].uv[uindex[2]]);
		if (t_dist < dist) {e1_index = cPoint->v2; e2_index = cPoint->v3; edge1_index = uindex[1]; edge2_index = uindex[2]; dist = t_dist;}

		// Dist to third edge
		t_dist = dist_to_line_segment_v2(pixel, tface[cPoint->index].uv[uindex[2]], tface[cPoint->index].uv[uindex[0]]);
		if (t_dist < dist) {e1_index = cPoint->v3; e2_index = cPoint->v1;  edge1_index = uindex[2]; edge2_index = uindex[0]; dist = t_dist;}


		/*
		*	Now find another face that is linked to that edge
		*/

		target_face = -1;

		for (i=0; i<numOfFaces; i++) {
			/*
			*	Check if both edge vertices share this face
			*/
			int v4 = -1;
			if (mface[i].v4) v4 = mface[i].v4;

			if ((e1_index == mface[i].v1 || e1_index == mface[i].v2 || e1_index == mface[i].v3 || e1_index == v4) &&
				(e2_index == mface[i].v1 || e2_index == mface[i].v2 || e2_index == mface[i].v3 || e2_index == v4)) {
				if (i == cPoint->index) continue;

				target_face = i;

				/*
				*	Get edge UV index
				*/
				if (e1_index == mface[i].v1) target_uv1 = 0;
				else if (e1_index == mface[i].v2) target_uv1 = 1;
				else if (e1_index == mface[i].v3) target_uv1 = 2;
				else target_uv1 = 3;

				if (e2_index == mface[i].v1) target_uv2 = 0;
				else if (e2_index == mface[i].v2) target_uv2 = 1;
				else if (e2_index == mface[i].v3) target_uv2 = 2;
				else target_uv2 = 3;

				break;
			}
		}

		//printf("target face loop done\n");
		// If none found return -1
		if (target_face == -1) return -1;

		/*
		*	If target face is connected in UV space as well, just use original index
		*/
		s_uv1 = (float *)tface[cPoint->index].uv[edge1_index];
		s_uv2 = (float *)tface[cPoint->index].uv[edge2_index];
		t_uv1 = (float *)tface[target_face].uv[target_uv1];
		t_uv2 = (float *)tface[target_face].uv[target_uv2];

		//printf("connected UV : %f,%f & %f,%f - %f,%f & %f,%f\n", s_uv1[0], s_uv1[1], s_uv2[0], s_uv2[1], t_uv1[0], t_uv1[1], t_uv2[0], t_uv2[1]);
		//if (s_uv1[0] == t_uv1[0]) printf("match 1\n");
		//if (s_uv1[0] == t_uv1[0] && s_uv1[1] == t_uv1[1]) printf("match 2\n");

		if (((s_uv1[0] == t_uv1[0] && s_uv1[1] == t_uv1[1]) &&
			 (s_uv2[0] == t_uv2[0] && s_uv2[1] == t_uv2[1]) ) ||
			((s_uv2[0] == t_uv1[0] && s_uv2[1] == t_uv1[1]) &&
			 (s_uv1[0] == t_uv2[0] && s_uv1[1] == t_uv2[1]) )) return ((px+neighX[n_index]) + surface->w*(py+neighY[n_index]));

		/*
		*	Find a point that is relatively at same edge position
		*	on this other face UV
		*/

		lambda = closest_to_line_v2(closest_point, pixel, tface[cPoint->index].uv[edge1_index], tface[cPoint->index].uv[edge2_index]);
		if (lambda < 0.0f) lambda = 0.0f;
		if (lambda > 1.0f) lambda = 1.0f;

		sub_v2_v2v2(dir_vec, tface[target_face].uv[target_uv2], tface[target_face].uv[target_uv1]);
		//normalize_v2(dir_vec);

		mul_v2_fl(dir_vec, lambda);

		copy_v2_v2(pixel, tface[target_face].uv[target_uv1]);
		add_v2_v2(pixel, dir_vec);
		pixel[0] = (pixel[0] * (float)surface->w) - 0.5f;
		pixel[1] = (pixel[1] * (float)surface->h) - 0.5f;

		final_pixel[0] = (int)floor(pixel[0]);
		final_pixel[1] = (int)floor(pixel[1]);

		//printf("target face %i, final pixel result %i, %i (%f, %f)\n", target_face, final_pixel[0], final_pixel[1], pixel[0], pixel[1]);

		if (final_pixel[0] < 0 || final_pixel[0] >= surface->w) return -1;
		if (final_pixel[1] < 0 || final_pixel[1] >= surface->h) return -1;

		final_index = final_pixel[0] + surface->w * final_pixel[1];

		if (surface->point[final_index].index != target_face) {/*printf("neigh point found index was not same!!! (%i, %i) -> (%i, %i) ... %i, %i\nedges : (%i, %i) -> (%i, %i) and\n(%i, %i) -> (%i, %i)\n", px+neighX[n_index], surface->h-1-(py+neighY[n_index]), final_pixel[0], surface->h-1-final_pixel[1], e1_index, target_uv1, 
			(int)floor(tface[cPoint->index].uv[edge1_index][0]*surface->w),surface->h-1-(int)floor(tface[cPoint->index].uv[edge1_index][1]*surface->h),(int)floor(tface[cPoint->index].uv[edge2_index][0]*surface->w),surface->h-1-(int)floor(tface[cPoint->index].uv[edge2_index][1]*surface->h),
			(int)floor(tface[target_face].uv[target_uv1][0]*surface->w),surface->h-1-(int)floor(tface[target_face].uv[target_uv1][1]*surface->h),(int)floor(tface[target_face].uv[target_uv2][0]*surface->w),surface->h-1-(int)floor(tface[target_face].uv[target_uv2][1]*surface->h));*/
			return -1;}

		//printf("final index found of %i\n", final_index);

		/*
		*	If final point is an "edge pixel", use it's "real" neighbour instead
		*/
		if (surface->point[final_index].neighbour_pixel != -1) final_index = cPoint->neighbour_pixel;

		return final_index;
		//closest_to_line_segment_v2(closest_point, float p[2], float l1[2], float l2[2])
	}
}

/*
*	Create Canvas Surface for baking
*/
int dynamicPaint_createCanvasSurface(DynamicPaintCanvasSettings *canvas) {

		int yy;
		int w,h;

		// Antialias jitter point relative coords
		float jitter5sample[10] =  {0.0f, 0.0f,
								-0.2f, -0.4f,
								0.2f, 0.4f,
								0.4f, -0.2f,
								-0.4f, 0.3f};


		struct DerivedMesh *dm = canvas->dm;
		int numOfFaces;
		MVert *mvert = NULL;
		MFace *mface = NULL;
		MTFace *tface = NULL;

		struct PaintSurface *surface;
		struct BB2d *faceBB = NULL;


		if (!dm) {printf("DynamicPaint bake failed: Invalid canvas mesh.\n"); return 0;}
		numOfFaces = dm->getNumFaces(dm);

		// Allocate memory for surface
		canvas->surface = (struct PaintSurface *) MEM_mapallocN(sizeof(struct PaintSurface), "MPCanvasSurface");
		surface = canvas->surface;
		surface->point = NULL;

		mvert = dm->getVertArray(dm);
		mface = dm->getFaceArray(dm);
		tface = DM_get_face_data_layer(dm, CD_MTFACE);

		printf("number of verts %i\n", dm->getNumVerts(dm));

		// Check for valid values
		if (!tface) {printf("DynamicPaint bake failed: No UV coordinates on canvas.\n");return 0;}
		if (canvas->resolution < 16 || canvas->resolution > 8096) {printf("DynamicPaint bake failed: Invalid texture resolution.\n"); return 0;}
	
		w = h = canvas->resolution;
		surface->w = w;
		surface->h = h;


		/*
		*	If effects enabled, calculate per face neighbour indexes
		*/
		if (canvas->effect) {


		}


		/*
		*	Start generating the surface
		*/
		printf("DynamicPaint: Preparing canvas of %ix%i pixels and %i faces.\n", w, h, numOfFaces);

		surface->pixelSamples = (canvas->flags & MOD_DPAINT_ANTIALIAS) ? 5 : 1;
		surface->point = (struct PaintSurfacePoint *) MEM_mapallocN(w*h*sizeof(struct PaintSurfacePoint), "PaintSurfaceData");


		/*
		*	Calculate face UV bonding boxes to optimize the search
		*/
		faceBB = (struct BB2d *) MEM_mapallocN(numOfFaces*sizeof(struct BB2d), "MPCanvasFaceBB");

		for (yy=0; yy<numOfFaces; yy++) {

			int numOfVert = (mface[yy].v4) ? 4 : 3;
			int i;
	
			VECCOPY2D(faceBB[yy].min, tface[yy].uv[0]);
			VECCOPY2D(faceBB[yy].max, tface[yy].uv[0]);

			for (i = 1; i<numOfVert; i++) {
				
				if (tface[yy].uv[i][0] < faceBB[yy].min[0]) faceBB[yy].min[0] = tface[yy].uv[i][0];
				if (tface[yy].uv[i][1] < faceBB[yy].min[1]) faceBB[yy].min[1] = tface[yy].uv[i][1];

				if (tface[yy].uv[i][0] > faceBB[yy].max[0]) faceBB[yy].max[0] = tface[yy].uv[i][0];
				if (tface[yy].uv[i][1] > faceBB[yy].max[1]) faceBB[yy].max[1] = tface[yy].uv[i][1];

			}
		}	// end face loop

		/*
		*	Allocate antialias sample data (without threads, coz malloc)
		*	(Non threadable)
		*/
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;
				PaintSurfacePoint *cPoint = (&surface->point[index]);

				// Initialize barycentricWeights
				cPoint->barycentricWeights = (struct Vec3f *) malloc( surface->pixelSamples * sizeof(struct Vec3f ));

			}
		} // end pixel loop



		/*
		*	Loop through every pixel and check
		*	if it's on domain mesh.
		*/
		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int i, sample;
				int index = xx+w*yy;
				PaintSurfacePoint *cPoint = (&surface->point[index]);

				short isInside = 0;	// if point is inside a face

				float d1[2], d2[2], d3[2], point[5][2];
				float dot00,dot01,dot02,dot11,dot12, invDenom, u,v;

				/*
				*	Init per pixel settings
				*/
				cPoint->color[0] = 0.0f;
				cPoint->color[1] = 0.0f;
				cPoint->color[2] = 0.0f;
				cPoint->alpha = 0.0f;
				cPoint->depth = 0.0f;

				cPoint->wetness = 0.0f;
				cPoint->e_alpha = 0.0f;
				cPoint->e_color[0] = 0.0f;
				cPoint->e_color[1] = 0.0f;
				cPoint->e_color[2] = 0.0f;
				cPoint->state = 0;

				cPoint->index = -1;
				cPoint->neighbour_pixel = -1;

				/*
				* A pixel middle sample isn't enough to find very narrow polygons
				* So using 4 samples on each corner too
				*/
				point[0][0] = ((float)xx) / w;
				point[0][1] = ((float)yy) / h;

				point[1][0] = ((float)xx+1) / w;
				point[1][1] = ((float)yy) / h;

				point[2][0] = ((float)xx) / w;
				point[2][1] = ((float)yy+1) / h;

				point[3][0] = ((float)xx+1) / w;
				point[3][1] = ((float)yy+1) / h;

				// Actual pixel center, used when collision is found
				point[4][0] = ((float)xx + 0.5f) / w;
				point[4][1] = ((float)yy + 0.5f) / h;


				// Loop through every face in the mesh
				for (i=0; i<numOfFaces; i++) {

						// Check uv bb
						if (faceBB[i].min[0] > (point[3][0])) continue;
						if (faceBB[i].min[1] > (point[3][1])) continue;
						if (faceBB[i].max[0] < (point[0][0])) continue;
						if (faceBB[i].max[1] < (point[0][1])) continue;
						
						for (sample=0; sample<5; sample++) {

							/*
							*	Calculate point inside a triangle check
							*	for uv0,1,2
							*/
							VECSUB2D(d1,  tface[i].uv[2], tface[i].uv[0]);	// uv2 - uv0
							VECSUB2D(d2,  tface[i].uv[1], tface[i].uv[0]);	// uv1 - uv0
							VECSUB2D(d3,  point[sample], tface[i].uv[0]);	// point - uv0


							dot00 = d1[0]*d1[0] + d1[1]*d1[1];
							dot01 = d1[0]*d2[0] + d1[1]*d2[1];
							dot02 = d1[0]*d3[0] + d1[1]*d3[1];
							dot11 = d2[0]*d2[0] + d2[1]*d2[1];
							dot12 = d2[0]*d3[0] + d2[1]*d3[1];

							invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
							u = (dot11 * dot02 - dot01 * dot12) * invDenom;
							v = (dot00 * dot12 - dot01 * dot02) * invDenom;

							if ((u > 0) && (v > 0) && (u + v < 1)) {isInside=1;} // is inside a triangle

							/*
							*	If collision wasn't found but the face is a quad
							*	do another check for the second half
							*/
							if ((!isInside) && mface[i].v4)
							{

								// change d2 to test the other half
								VECSUB2D(d2,  tface[i].uv[3], tface[i].uv[0]);	// uv3 - uv0

								// test again
								dot00 = d1[0]*d1[0] + d1[1]*d1[1];
								dot01 = d1[0]*d2[0] + d1[1]*d2[1];
								dot02 = d1[0]*d3[0] + d1[1]*d3[1];
								dot11 = d2[0]*d2[0] + d2[1]*d2[1];
								dot12 = d2[0]*d3[0] + d2[1]*d3[1];

								invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
								u = (dot11 * dot02 - dot01 * dot12) * invDenom;
								v = (dot00 * dot12 - dot01 * dot02) * invDenom;

								if ((u > 0) && (v > 0) && (u + v < 1)) {isInside=2;} // is inside a second half of quad

							}


							/*
							*	If point was inside the face
							*/
							if (isInside != 0) {

								float uv1co[2], uv2co[2], uv3co[2], uv[2];
								int j;

								// Get triagnle uvs
								if (isInside==1) {
									VECCOPY2D(uv1co, tface[i].uv[0]);
									VECCOPY2D(uv2co, tface[i].uv[1]);
									VECCOPY2D(uv3co, tface[i].uv[2]);
								}
								else {
									VECCOPY2D(uv1co, tface[i].uv[0]);
									VECCOPY2D(uv2co, tface[i].uv[2]);
									VECCOPY2D(uv3co, tface[i].uv[3]);
								}

								// Add b-weights per anti-aliasing sample
								for (j=0; j<surface->pixelSamples; j++) {

									uv[0] = point[4][0] + jitter5sample[j*2] / w;
									uv[1] = point[4][1] + jitter5sample[j*2+1] / h;

									barycentric_weights_v2(uv1co, uv2co, uv3co, uv, cPoint->barycentricWeights[j].v);
								}

								// Set surface point values for collision
								cPoint->index = i;							// face index
								cPoint->quad = (isInside == 2) ? 1 : 0;		// quad or tri

								// save vertex indexes
								cPoint->v1 = (isInside == 2) ? mface[i].v1 : mface[i].v1;
								cPoint->v2 = (isInside == 2) ? mface[i].v3 : mface[i].v2;
								cPoint->v3 = (isInside == 2) ? mface[i].v4 : mface[i].v3;
								
								i = numOfFaces;	// make sure we exit face loop as well
								break;
							}	// end isInside
						} // end sample;


				}	// end face loop

			}	// end of yy loop
		}	// end of xx loop




		/*
		*	Now loop through every pixel that was left without index
		*	and find if they have neighbour pixels that have an index.
		*	If so use that polygon as pixel surface.
		*	(To avoid seams on uv island edges.)
		*/

		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;
				PaintSurfacePoint *cPoint = (&surface->point[index]);

				// If point isnt't on canvas mesh
				if (cPoint->index == -1) {
					int u_min, u_max, v_min, v_max;
					int u,v, ind;
					float point[2];

					// get loop area
					u_min = (xx > 0) ? -1 : 0;
					u_max = (xx < (w-1)) ? 1 : 0;
					v_min = (yy > 0) ? -1 : 0;
					v_max = (yy < (h-1)) ? 1 : 0;

					point[0] = ((float)xx + 0.5f) / w;
					point[1] = ((float)yy + 0.5f) / h;

					// search through defined area for neighbour
					for (u=u_min; u<=u_max; u++)
						for (v=v_min; v<=v_max; v++) {

							// if not this pixel itself
							if (u!=0 || v!=0) {
								ind = (xx+u)+w*(yy+v);

								// if neighbour has index
								if (surface->point[ind].index != -1) {

									float uv1co[2], uv2co[2], uv3co[2], uv[2];
									int i = surface->point[ind].index, j;

									/*
									*	Now calculate pixel data for this pixel as it was on polygon surface
									*/

									if (!surface->point[ind].quad) {
										VECCOPY2D(uv1co, tface[i].uv[0]);
										VECCOPY2D(uv2co, tface[i].uv[1]);
										VECCOPY2D(uv3co, tface[i].uv[2]);
									}
									else {
										VECCOPY2D(uv1co, tface[i].uv[0]);
										VECCOPY2D(uv2co, tface[i].uv[2]);
										VECCOPY2D(uv3co, tface[i].uv[3]);
									}

									// Add b-weights per anti-aliasing sample
									for (j=0; j<surface->pixelSamples; j++) {

										uv[0] = point[0] + jitter5sample[j*2] / w;
										uv[1] = point[1] + jitter5sample[j*2+1] / h;

										barycentric_weights_v2(uv1co, uv2co, uv3co, uv, cPoint->barycentricWeights[j].v);
									}

									// Set values
									cPoint->neighbour_pixel = ind;				// face index
									cPoint->quad = surface->point[ind].quad;		// quad or tri

									// save vertex indexes
									cPoint->v1 = (cPoint->quad) ? mface[i].v1 : mface[i].v1;
									cPoint->v2 = (cPoint->quad) ? mface[i].v3 : mface[i].v2;
									cPoint->v3 = (cPoint->quad) ? mface[i].v4 : mface[i].v3;

									u = u_max + 1;	// make sure we exit outer loop as well
									break;
								}

							} // end itself check
						} // end uv loop

				}	// end if has index

			} } // end pixel loop


		/*
		*	When base loop is over convert found neighbour indexes to real ones
		*	Also count the final number of active surface points
		*/

		surface->active_points = 0;

		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;
				PaintSurfacePoint *cPoint = (&surface->point[index]);

				if (cPoint->index == -1 && cPoint->neighbour_pixel != -1) cPoint->index = surface->point[cPoint->neighbour_pixel].index;

				if (cPoint->index != -1) surface->active_points++;

			}
		} // end pixel loop

#if 0
		/*
		*	-----------------------------------------------------------------
		*	For debug, output pixel statuses to color map
		*	-----------------------------------------------------------------
		*/

		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;
				PaintSurfacePoint *cPoint = (&surface->point[index]);

				// Every pixel that is assigned as "edge pixel" gets red color
				if (cPoint->neighbour_pixel != -1) {cPoint->color[2] = 1.0f; cPoint->alpha=1.0f;}
				// and every directly inside polygon aligned pixel gets blue color
				if (cPoint->index != -1) {cPoint->color[0] = 1.0f; cPoint->alpha=1.0f;}

			}
		} // end pixel loop

#endif



		/*
		*	If effects enabled, create surface effect / wet layer
		*	neighbour lists. Processes possibly moving data.
		*/
		if (canvas->effect) {

			#pragma omp parallel for schedule(static)
			for (yy = 0; yy < h; yy++)
			{
				int xx;
				for (xx = 0; xx < w; xx++)
				{
					int i;
					PaintSurfacePoint *cPoint = (&surface->point[xx+w*yy]);

					for (i=0; i<8; i++) {
						int x,y, index;
						PaintSurfacePoint *tPoint;

						cPoint->neighbours[i] = -1;

						if (cPoint->index != -1) {

							x = xx + neighX[i];
							y = yy + neighY[i];

							if (x<0 || x>=w) continue;
							if (y<0 || y>=h) continue;

							index = x+w*y;
							tPoint = (&surface->point[index]);

							// if target point is on same face, mark it as neighbour
							// (and if it isn't an "edge pixel")
							if ((tPoint->index == cPoint->index) && (tPoint->neighbour_pixel == -1)) {
								cPoint->neighbours[i] = index;
							}
							else {	// else find the "real" neighbour
								int neigh = dynamicPaint_findNeighbourPixel(canvas, xx, yy, i);

								if (neigh != -1) cPoint->neighbours[i] = neigh;
							}
						}
					}

				}
			} // end pixel loop
		} // effect


		MEM_freeN(faceBB);

		return 1;
}


/*
*	Free canvas surface
*/
void dynamicPaint_cleanCanvasSurface(DynamicPaintCanvasSettings *canvas) {

	int w,h,k;
	
	if (!canvas) return;
	if (!canvas->surface) return;

	w = canvas->surface->w;
	h = canvas->surface->h;

	#pragma omp parallel for schedule(static,1)
	for (k = 0; k < w*h; k++)
	{
		free(canvas->surface->point[k].barycentricWeights);
	}

	if (canvas->surface->point) MEM_freeN(canvas->surface->point);
	MEM_freeN(canvas->surface);
}

/* A modified callback to bvh tree raycast. The tree must bust have been built using bvhtree_from_mesh_faces.
/* userdata must be a BVHMeshCallbackUserdata built from the same mesh as the tree.
*  
*	To optimize paint detection speed this doesn't calculate hit coordinates or normal.
*	If ray hit the second half of a quad, no[0] is set to 1.0f.
*/
static void mesh_faces_spherecast_dp(void *userdata, int index, const BVHTreeRay *ray, BVHTreeRayHit *hit)
{
	const BVHTreeFromMesh *data = (BVHTreeFromMesh*) userdata;
	MVert *vert	= data->vert;
	MFace *face = data->face + index;
	short quad = 0;

	float *t0, *t1, *t2, *t3;
	t0 = vert[ face->v1 ].co;
	t1 = vert[ face->v2 ].co;
	t2 = vert[ face->v3 ].co;
	t3 = face->v4 ? vert[ face->v4].co : NULL;

	
	do
	{	
		float dist;
		dist = ray_tri_intersection(ray, hit->dist, t0, t1, t2);

		if(dist >= 0 && dist < hit->dist)
		{
			hit->index = index;
			hit->dist = dist;
			//VECADDFAC(hit->co, ray->origin, ray->direction, dist);

			//normal_tri_v3( hit->no,t0, t1, t2);
			hit->no[0] = (quad) ? 1.0f : 0.0f;
		}

		t1 = t2;
		t2 = t3;
		t3 = NULL;
		quad = 1;

	} while(t2);
}

/* A modified callback to bvh tree nearest point. The tree must bust have been built using bvhtree_from_mesh_faces.
*  userdata must be a BVHMeshCallbackUserdata built from the same mesh as the tree.
*  
*	To optimize paint detection speed this doesn't calculate hit normal.
*	If ray hit the second half of a quad, no[0] is set to 1.0f, else 0.0f
*/
static void mesh_faces_nearest_point_dp(void *userdata, int index, const float *co, BVHTreeNearest *nearest)
{
	const BVHTreeFromMesh *data = (BVHTreeFromMesh*) userdata;
	MVert *vert	= data->vert;
	MFace *face = data->face + index;
	short quad = 0;

	float *t0, *t1, *t2, *t3;
	t0 = vert[ face->v1 ].co;
	t1 = vert[ face->v2 ].co;
	t2 = vert[ face->v3 ].co;
	t3 = face->v4 ? vert[ face->v4].co : NULL;

	
	do
	{	
		float nearest_tmp[3], dist;
		int vertex, edge;
		
		dist = nearest_point_in_tri_surface(t0, t1, t2, co, &vertex, &edge, nearest_tmp);
		if(dist < nearest->dist)
		{
			nearest->index = index;
			nearest->dist = dist;
			VECCOPY(nearest->co, nearest_tmp);
			//normal_tri_v3( nearest->no,t0, t1, t2);
			nearest->no[0] = (quad) ? 1.0f : 0.0f;
		}

		t1 = t2;
		t2 = t3;
		t3 = NULL;
		quad = 1;

	} while(t2);
}


/*
*	Calculate inverse matrices for material related objects
*	in case texture mapping is object related
*	(obj->imat isn't auto-updated)
*/
void DynamicPaint_InitMaterialObjects(Object *paintOb) {

	Material *mat = give_current_material(paintOb, paintOb->actcol);
	MTex *mtex;
	Tex *tex;
	int tex_nr;

	/*
	*	Calculate inverse transformation matrix
	*	for this object
	*/
	invert_m4_m4(paintOb->imat, paintOb->obmat);

	if (mat == NULL) {printf("mat was null...\n");return;}

	/*
	*	Loop through every material texture and check
	*	if they are mapped by other object
	*/
	
	for(tex_nr=0; tex_nr<MAX_MTEX; tex_nr++) {
		/* separate tex switching */
		if(mat->septex & (1<<tex_nr)) continue;
	
		if(mat->mtex[tex_nr]) {
			mtex= mat->mtex[tex_nr];
			tex= mtex->tex;

			if(tex==0) continue;
			
			/* which coords */
			if(mtex->texco==TEXCO_OBJECT) { 
				Object *ob= mtex->object;
				if(ob) {						
					invert_m4_m4(ob->imat, ob->obmat);
				}
			}

		}
	}	// end texture loop

}

// USE MFace mat_nr to define material!!!		-----------------

#if 1

/*
*	Edited version of do_material_tex()

*	also see shade_input_set_shade_texco() for ORCO settings
*		 and shade_input_set_uv() for face u,v calculation
*/

/* a modified part of shadeinput.c -> shade_input_set_uv() / shade_input_set_shade_texco() */
void textured_face_generate_uv(float *uv, float *normal, float *hit, float *v1, float *v2, float *v3)
{

	/* most of this could become re-used for faces */
	float detsh, t00, t10, t01, t11, xn, yn, zn;
	int axis1, axis2;

	/* find most stable axis to project */
	xn= fabs(normal[0]);
	yn= fabs(normal[1]);
	zn= fabs(normal[2]);

	if(zn>=xn && zn>=yn) { axis1= 0; axis2= 1; }
	else if(yn>=xn && yn>=zn) { axis1= 0; axis2= 2; }
	else { axis1= 1; axis2= 2; }

	/* compute u,v and derivatives */
	t00= v3[axis1]-v1[axis1]; t01= v3[axis2]-v1[axis2];
	t10= v3[axis1]-v2[axis1]; t11= v3[axis2]-v2[axis2];

	detsh= 1.0f/(t00*t11-t10*t01);
	t00*= detsh; t01*=detsh; 
	t10*=detsh; t11*=detsh;

	uv[0] = (hit[axis1]-v3[axis1])*t11-(hit[axis2]-v3[axis2])*t10;
	uv[1] = (hit[axis2]-v3[axis2])*t00-(hit[axis1]-v3[axis1])*t01;

	/* u and v are in range -1 to 0, we allow a little bit extra but not too much, screws up speedvectors */
	CLAMP(uv[0], -2.0f, 1.0f);
	CLAMP(uv[1], -2.0f, 1.0f);
}

/* a modified part of shadeinput.c -> shade_input_set_uv() / shade_input_set_shade_texco() */
void textured_face_get_uv(float *uv_co, float *normal, float *uv, int faceIndex, short quad, MTFace *tface)
{
	float *uv1, *uv2, *uv3;
	float l;

	l= 1.0f+uv[0]+uv[1];
		
	uv1= tface[faceIndex].uv[0];
	uv2= (quad) ? tface[faceIndex].uv[2] : tface[faceIndex].uv[1];
	uv3= (quad) ? tface[faceIndex].uv[3] : tface[faceIndex].uv[2];
				
	uv_co[0]= -1.0f + 2.0f*(l*uv3[0]-uv[0]*uv1[0]-uv[1]*uv2[0]);
	uv_co[1]= -1.0f + 2.0f*(l*uv3[1]-uv[0]*uv1[1]-uv[1]*uv2[1]);
	uv_co[2]= 0.0f;	/* texture.c assumes there are 3 coords */
}

void DynamicPaint_SampleSolidMaterial(float color[3], float *alpha, Material *mat, Object *paintOb, float xyz[3], int faceIndex, short isQuad, DerivedMesh *orcoDm)
{
	MTex *mtex;
	Tex *tex;
	TexResult texres= {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, NULL};
	float co[3], xyz_local[3];
	float fact, stencilTin=1.0;
	float texvec[3];
	int tex_nr, rgbnor= 0, warpdone=0;

	float uv[3], normal[3];

	MFace *mface;
	int v1, v2, v3;
	MVert *mvert;
	
	// Get face data

	mvert = orcoDm->getVertArray(orcoDm);
	mface = orcoDm->getFaceArray(orcoDm);
	v1=mface[faceIndex].v1, v2=mface[faceIndex].v2, v3=mface[faceIndex].v3;
	if (isQuad) {v2=mface[faceIndex].v3; v3=mface[faceIndex].v4;}

	normal_tri_v3( normal, mvert[v1].co, mvert[v2].co, mvert[v3].co);

	// Assign material base values
	color[0] = mat->r;
	color[1] = mat->g;
	color[2] = mat->b;
	*alpha = mat->alpha;

	VECCOPY(xyz_local, xyz);
	mul_m4_v3(paintOb->imat, xyz_local);

	for(tex_nr=0; tex_nr<MAX_MTEX; tex_nr++) {
		
		/* separate tex switching */
			if(mat->septex & (1<<tex_nr)) continue;
		
			if(mat->mtex[tex_nr]) {
				mtex= mat->mtex[tex_nr];
				tex= mtex->tex;
			
			tex= mtex->tex;
			if(tex==0) continue;

			VECCOPY(co, xyz);

			/* which coords */
			if(mtex->texco==TEXCO_ORCO) {
					float l;

					/*
					*	Get generated UV
					*/
					textured_face_generate_uv(uv, normal, xyz_local, mvert[v1].co, mvert[v2].co, mvert[v3].co);

					l= 1.0f+uv[0]+uv[1];

					// calculate generated coordinate
					// ** Keep up-to-date with shadeinput.c -> shade_input_set_shade_texco() **
					co[0]= l*mvert[v3].co[0]-uv[0]*mvert[v1].co[0]-uv[1]*mvert[v2].co[0];
					co[1]= l*mvert[v3].co[1]-uv[0]*mvert[v1].co[1]-uv[1]*mvert[v2].co[1];
					co[2]= l*mvert[v3].co[2]-uv[0]*mvert[v1].co[2]-uv[1]*mvert[v2].co[2];
					// 
					// shade_input_set_shade_texco
			}
			// Ignore Sticky
			else if(mtex->texco==TEXCO_OBJECT) {
				Object *ob= mtex->object;

				VECCOPY(co, xyz);
				// convert from world space to paint space
				mul_m4_v3(paintOb->imat, co);
				if(ob) {
					mul_m4_v3(ob->imat, co);
				}
			}
			else if(mtex->texco==TEXCO_REFL) {
				//co= shi->ref; dx= shi->dxref; dy= shi->dyref;
			}
			else if(mtex->texco==TEXCO_NORM) {
				//co= shi->orn; dx= shi->dxno; dy= shi->dyno;
			}
			else if(mtex->texco==TEXCO_TANGENT) {
				//co= shi->tang; dx= shi->dxno; dy= shi->dyno;
			}
			else if(mtex->texco==TEXCO_GLOB) {
				VECCOPY(co, xyz);
			}
			else if(mtex->texco==TEXCO_UV) {
				MTFace *tface;

				if(mtex->uvname[0] != 0) {
					/*for(i = 0; i < shi->totuv; i++) {
						if(strcmp(shi->uv[i].name, mtex->uvname)==0) {
							suv= &shi->uv[i];
							break;
						}
					}*/
					tface = CustomData_get_layer_named(&orcoDm->faceData, CD_MTFACE, mtex->uvname);
				}
				else tface = DM_get_face_data_layer(orcoDm, CD_MTFACE);

				// Get generated coordinates to calculate UV from
				textured_face_generate_uv(uv, normal, xyz_local, mvert[v1].co, mvert[v2].co, mvert[v3].co);

				/*
				*	Get UV mapping coordinate
				*/
				textured_face_get_uv(co, normal, uv, faceIndex, isQuad, tface);
			}
			else if(mtex->texco==TEXCO_WINDOW) {
				//co= shi->winco; dx= shi->dxwin; dy= shi->dywin;
			}
			else if(mtex->texco==TEXCO_STRAND) {
				//
			}
			else if(mtex->texco==TEXCO_STRESS) {
				//
			}
			else continue;	// can happen when texco defines disappear and it renders old files

			{
				texco_mapping_ext(normal, tex, mtex, co, 0, 0, texvec);

				if(tex->use_nodes && tex->nodetree) {
					// No support for nodes yet.
					continue;
				}
				else
					//rgbnor = multitex(mtex->tex, co, 0, 0, 0, &texres, 0, 0);
					rgbnor = multitex_ext(mtex->tex, co, 0, 0, 0, &texres);

			}

			/* texture output */

			if( (rgbnor & TEX_RGB) && (mtex->texflag & MTEX_RGBTOINT)) {
				texres.tin= (0.35*texres.tr+0.45*texres.tg+0.2*texres.tb);
				rgbnor-= TEX_RGB;
			}
			if(mtex->texflag & MTEX_NEGATIVE) {
				if(rgbnor & TEX_RGB) {
					texres.tr= 1.0-texres.tr;
					texres.tg= 1.0-texres.tg;
					texres.tb= 1.0-texres.tb;
				}
				texres.tin= 1.0-texres.tin;
			}
			if(mtex->texflag & MTEX_STENCIL) {
				if(rgbnor & TEX_RGB) {
					fact= texres.ta;
					texres.ta*= stencilTin;
					stencilTin*= fact;
				}
				else {
					fact= texres.tin;
					texres.tin*= stencilTin;
					stencilTin*= fact;
				}
			}
			

			/* mapping */
			if(mtex->mapto & (MAP_COL)) {
				float tcol[3];
				
				/* stencil maps on the texture control slider, not texture intensity value */
				
				tcol[0]=texres.tr; tcol[1]=texres.tg; tcol[2]=texres.tb;
				
				if((rgbnor & TEX_RGB)==0) {
					tcol[0]= mtex->r;
					tcol[1]= mtex->g;
					tcol[2]= mtex->b;
				}
				else if(mtex->mapto & MAP_ALPHA) {
					texres.tin= stencilTin;
				}
				else texres.tin= texres.ta;
				
				/* inverse gamma correction */
				if (tex->type==TEX_IMAGE) {
					Image *ima = tex->ima;
					ImBuf *ibuf = BKE_image_get_ibuf(ima, &tex->iuser);
					
					/* don't linearize float buffers, assumed to be linear */
					/*if (ibuf && !(ibuf->rect_float) && R.r.color_mgt_flag & R_COLOR_MANAGEMENT)
						srgb_to_linearrgb_v3_v3(tcol, tcol);*/
				}
				
				if(mtex->mapto & MAP_COL) {
					float colfac= mtex->colfac*stencilTin;
					texture_rgb_blend(color, tcol, color, texres.tin, colfac, mtex->blendtype);
				}
				/*if(mtex->mapto & MAP_COLSPEC) {
					float colspecfac= mtex->colspecfac*stencilTin;
					texture_rgb_blend(&shi->specr, tcol, &shi->specr, texres.tin, colspecfac, mtex->blendtype);
				}
				if(mtex->mapto & MAP_COLMIR) {
					float mirrfac= mtex->mirrfac*stencilTin;

					// exception for envmap only
					if(tex->type==TEX_ENVMAP && mtex->blendtype==MTEX_BLEND) {
						fact= texres.tin*mirrfac;
						facm= 1.0- fact;
						shi->refcol[0]= fact + facm*shi->refcol[0];
						shi->refcol[1]= fact*tcol[0] + facm*shi->refcol[1];
						shi->refcol[2]= fact*tcol[1] + facm*shi->refcol[2];
						shi->refcol[3]= fact*tcol[2] + facm*shi->refcol[3];
					}
					else {
						texture_rgb_blend(&shi->mirr, tcol, &shi->mirr, texres.tin, mirrfac, mtex->blendtype);
					}
				}*/
			}
			

			if(mtex->mapto & MAP_VARS) {
				/* stencil maps on the texture control slider, not texture intensity value */
				
				if(rgbnor & TEX_RGB) {
					if(texres.talpha) texres.tin= texres.ta;
					else texres.tin= (0.35*texres.tr+0.45*texres.tg+0.2*texres.tb);
				}

				/*if(mtex->mapto & MAP_REF) {
					float difffac= mtex->difffac*stencilTin;

					shi->refl= texture_value_blend(mtex->def_var, shi->refl, texres.tin, difffac, mtex->blendtype);
					if(shi->refl<0.0) shi->refl= 0.0;
				}
				if(mtex->mapto & MAP_SPEC) {
					float specfac= mtex->specfac*stencilTin;
					
					shi->spec= texture_value_blend(mtex->def_var, shi->spec, texres.tin, specfac, mtex->blendtype);
					if(shi->spec<0.0) shi->spec= 0.0;
				}
				if(mtex->mapto & MAP_EMIT) {
					float emitfac= mtex->emitfac*stencilTin;

					shi->emit= texture_value_blend(mtex->def_var, shi->emit, texres.tin, emitfac, mtex->blendtype);
					if(shi->emit<0.0) shi->emit= 0.0;
				}*/
				if(mtex->mapto & MAP_ALPHA) {
					float alphafac= mtex->alphafac*stencilTin;

					*alpha= texture_value_blend(mtex->def_var, *alpha, texres.tin, alphafac, mtex->blendtype);
					if(*alpha<0.0) *alpha= 0.0;
					else if(*alpha>1.0) *alpha= 1.0;
				}
				/*if(mtex->mapto & MAP_HAR) {
					float har;  // have to map to 0-1
					float hardfac= mtex->hardfac*stencilTin;
					
					har= ((float)shi->har)/128.0;
					har= 128.0*texture_value_blend(mtex->def_var, har, texres.tin, hardfac, mtex->blendtype);
					
					if(har<1.0) shi->har= 1; 
					else if(har>511.0) shi->har= 511;
					else shi->har= (int)har;
				}
				if(mtex->mapto & MAP_RAYMIRR) {
					float raymirrfac= mtex->raymirrfac*stencilTin;

					shi->ray_mirror= texture_value_blend(mtex->def_var, shi->ray_mirror, texres.tin, raymirrfac, mtex->blendtype);
					if(shi->ray_mirror<0.0) shi->ray_mirror= 0.0;
					else if(shi->ray_mirror>1.0) shi->ray_mirror= 1.0;
				}
				if(mtex->mapto & MAP_TRANSLU) {
					float translfac= mtex->translfac*stencilTin;

					shi->translucency= texture_value_blend(mtex->def_var, shi->translucency, texres.tin, translfac, mtex->blendtype);
					if(shi->translucency<0.0) shi->translucency= 0.0;
					else if(shi->translucency>1.0) shi->translucency= 1.0;
				}
				if(mtex->mapto & MAP_AMB) {
					float ambfac= mtex->ambfac*stencilTin;

					shi->amb= texture_value_blend(mtex->def_var, shi->amb, texres.tin, ambfac, mtex->blendtype);
					if(shi->amb<0.0) shi->amb= 0.0;
					else if(shi->amb>1.0) shi->amb= 1.0;
					
					shi->ambr= shi->amb*R.wrld.ambr;
					shi->ambg= shi->amb*R.wrld.ambg;
					shi->ambb= shi->amb*R.wrld.ambb;
				}*/
			}
		}
	}
}
#endif

void DynamicPaint_SampleVolumeMaterial(float color[3], float *alpha, Material *mat, Object *paintOb, float xyz[3]) {
		int mapto_flag  = MAP_DENSITY | MAP_REFLECTION_COL;
		float *col = color;

		color[0] = mat->vol.reflection_col[0];
		color[1] = mat->vol.reflection_col[1];
		color[2] = mat->vol.reflection_col[2];
		*alpha = mat->vol.density;

		//printf("material color: %f, %f, %f\n",mat->r, mat->g, mat->b);
		//printf("material color direct: %f, %f, %f\n",paintOb->mat[paintOb->actcol].r, paintOb->mat[paintOb->actcol].g, paintOb->mat[paintOb->actcol].b);

		if (1)
		{
			MTex *mtex;
			Tex *tex;
			TexResult texres= {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, NULL};
			int tex_nr, rgbnor= 0;
			float co[3], texvec[3];
			float fact, stencilTin=1.0;
	
			for(tex_nr=0; tex_nr<MAX_MTEX; tex_nr++) {
				/* separate tex switching */
				if(mat->septex & (1<<tex_nr)) continue;
		
				if(mat->mtex[tex_nr]) {
					mtex= mat->mtex[tex_nr];
					tex= mtex->tex;

					if(tex==0) continue;

					if (1) {

			
					/* only process if this texture is mapped 
					 * to one that we're interested in */
					if (!(mtex->mapto & mapto_flag)) continue;
			
					/* which coords */
					if(mtex->texco==TEXCO_OBJECT) { 
						Object *ob= mtex->object;
						ob= mtex->object;
						if(ob) {						
							VECCOPY(co, xyz);
							mul_m4_v3(ob->imat, co);
						}
					}
					/* not really orco, but 'local' */
					else if(mtex->texco==TEXCO_ORCO) {
				
						{
							Object *ob= paintOb;
							VECCOPY(co, xyz);
							mul_m4_v3(ob->imat, co);

							//printf("paint obj orig: %f, %f, %f\n", ob->orig[0], ob->orig[1], ob->orig[2]);
						}
					}
					else if(mtex->texco==TEXCO_GLOB) {							
					   VECCOPY(co, xyz);
					   //mul_m4_v3(R.viewinv, co);
					}
					else continue;	// can happen when texco defines disappear and it renders old files

					texres.nor= NULL;
			
					if(tex->type==TEX_IMAGE) {
						continue;	/* not supported yet */				
						//do_2d_mapping(mtex, texvec, NULL, NULL, dxt, dyt);
					}
					else {
						/* placement */
						if(mtex->projx) texvec[0]= mtex->size[0]*(co[mtex->projx-1]+mtex->ofs[0]);
						else texvec[0]= mtex->size[0]*(mtex->ofs[0]);

						if(mtex->projy) texvec[1]= mtex->size[1]*(co[mtex->projy-1]+mtex->ofs[1]);
						else texvec[1]= mtex->size[1]*(mtex->ofs[1]);

						if(mtex->projz) texvec[2]= mtex->size[2]*(co[mtex->projz-1]+mtex->ofs[2]);
						else texvec[2]= mtex->size[2]*(mtex->ofs[2]);
					}
					}
			
					rgbnor= multitex_ext(tex, texvec, NULL, NULL, 0, &texres);
			
					/* texture output */

					if( (rgbnor & TEX_RGB) && (mtex->texflag & MTEX_RGBTOINT)) {
						texres.tin= (0.35*texres.tr+0.45*texres.tg+0.2*texres.tb);
						rgbnor-= TEX_RGB;
					}
					if(mtex->texflag & MTEX_NEGATIVE) {
						if(rgbnor & TEX_RGB) {
							texres.tr= 1.0-texres.tr;
							texres.tg= 1.0-texres.tg;
							texres.tb= 1.0-texres.tb;
						}
						texres.tin= 1.0-texres.tin;
					}
					if(mtex->texflag & MTEX_STENCIL) {
						if(rgbnor & TEX_RGB) {
							fact= texres.ta;
							texres.ta*= stencilTin;
							stencilTin*= fact;
						}
						else {
							fact= texres.tin;
							texres.tin*= stencilTin;
							stencilTin*= fact;
						}
					}
			

					if((mapto_flag & (MAP_EMISSION_COL+MAP_TRANSMISSION_COL+MAP_REFLECTION_COL)) && (mtex->mapto & (MAP_EMISSION_COL+MAP_TRANSMISSION_COL+MAP_REFLECTION_COL))) {
						float tcol[3];
				
						/* stencil maps on the texture control slider, not texture intensity value */
				
						if((rgbnor & TEX_RGB)==0) {
							tcol[0]= mtex->r;
							tcol[1]= mtex->g;
							tcol[2]= mtex->b;
						} else {
							tcol[0]=texres.tr;
							tcol[1]=texres.tg;
							tcol[2]=texres.tb;
							if(texres.talpha)
								texres.tin= texres.ta;
						}
				
						/* used for emit */
						if((mapto_flag & MAP_EMISSION_COL) && (mtex->mapto & MAP_EMISSION_COL)) {
							float colemitfac= mtex->colemitfac*stencilTin;
							texture_rgb_blend(col, tcol, col, texres.tin, colemitfac, mtex->blendtype);
						}
				
						if((mapto_flag & MAP_REFLECTION_COL) && (mtex->mapto & MAP_REFLECTION_COL)) {
							float colreflfac= mtex->colreflfac*stencilTin;
							texture_rgb_blend(col, tcol, col, texres.tin, colreflfac, mtex->blendtype);
						}
				
						if((mapto_flag & MAP_TRANSMISSION_COL) && (mtex->mapto & MAP_TRANSMISSION_COL)) {
							float coltransfac= mtex->coltransfac*stencilTin;
							texture_rgb_blend(col, tcol, col, texres.tin, coltransfac, mtex->blendtype);
						}
					}
			
					if((mapto_flag & MAP_VARS) && (mtex->mapto & MAP_VARS)) {
						/* stencil maps on the texture control slider, not texture intensity value */
				
						/* convert RGB to intensity if intensity info isn't provided */
						if (!(rgbnor & TEX_INT)) {
							if (rgbnor & TEX_RGB) {
								if(texres.talpha) texres.tin= texres.ta;
								else texres.tin= (0.35*texres.tr+0.45*texres.tg+0.2*texres.tb);
							}
						}
						if((mapto_flag & MAP_DENSITY) && (mtex->mapto & MAP_DENSITY)) {
							float densfac= mtex->densfac*stencilTin;

							*alpha = texture_value_blend(mtex->def_var, *alpha, texres.tin, densfac, mtex->blendtype);
							CLAMP(*alpha, 0.0, 1.0);
						}
						/*if((mapto_flag & MAP_REFLECTION) && (mtex->mapto & MAP_REFLECTION)) {
							float reflfac= mtex->reflfac*stencilTin;
					
							*val = texture_value_blend(mtex->def_var, *val, texres.tin, reflfac, mtex->blendtype);
							CLAMP(*val, 0.0, 1.0);
						}*/
					}
				}
			}

		}
}

void DynamicPaint_GetMaterialColor(float *color, float *alpha, Object *paintOb, float pixelCoord[3], float paintHit[3], int faceIndex, short isQuad, DerivedMesh *orcoDm) {


	Material *material;
	MFace *mface;

	/*
	*	Get face material
	*/
	mface = orcoDm->getFaceArray(orcoDm);
	material = give_current_material(paintOb, mface[faceIndex].mat_nr+1);
	if (material == NULL) return;	// No material assigned


	/*
	*	Sample textured material color in given position depending on material type
	*/
	if (material->material_type == MA_TYPE_SURFACE) {

		DynamicPaint_SampleSolidMaterial(color, alpha, material, paintOb, paintHit, faceIndex, isQuad, orcoDm);
	}
	else if (material->material_type == MA_TYPE_VOLUME) {
		
		DynamicPaint_SampleVolumeMaterial(color, alpha, material, paintOb, pixelCoord);

	}
	else if (material->material_type == MA_TYPE_HALO) {
		
		// change nothing

	}
}

void DynamicPaint_MixPaintColors(PaintSurfacePoint *cPoint, int paintFlags, float *paintColor, float *paintAlpha, float *paintWetness)
{
	float invFact = 1.0f - (*paintAlpha);

	if (!(paintFlags & MOD_DPAINT_EREASE)) {

		// If this point has previous paint
		if (cPoint->e_alpha > 0)
		{
			/*
			*	Mix colors by the factor
			*/
			cPoint->e_color[0] = cPoint->e_color[0]*invFact + paintColor[0]*(*paintAlpha);
			cPoint->e_color[1] = cPoint->e_color[1]*invFact + paintColor[1]*(*paintAlpha);
			cPoint->e_color[2] = cPoint->e_color[2]*invFact + paintColor[2]*(*paintAlpha);
		}
		else
		{
			// set first color value straight to paint color
			cPoint->e_color[0] = paintColor[0];
			cPoint->e_color[1] = paintColor[1];
			cPoint->e_color[2] = paintColor[2];
		}

		// Increase alpha if higher than existing
		if (cPoint->e_alpha < (*paintAlpha)) cPoint->e_alpha = (*paintAlpha);
		if (cPoint->wetness < (*paintWetness)) cPoint->wetness = (*paintWetness);
	}
	else {
		float a_ratio, a_highest;
		float wetness = (1.0f - (*paintWetness)) * (*paintAlpha);

		/*
		*	Make highest alpha to match ereased value
		*	but maintain alpha ratio
		*/
		a_highest = (cPoint->e_alpha > cPoint->alpha) ? cPoint->e_alpha : cPoint->alpha;
		if (a_highest > invFact) {
			a_ratio = invFact / a_highest;

			cPoint->e_alpha *= a_ratio;
			cPoint->alpha *= a_ratio;
		}
		if (cPoint->wetness > wetness) cPoint->wetness = wetness;
	}
}

/*
*	Paint a mesh to canvas
*/
int DynamicPaint_PaintMesh(DynamicPaintCanvasSettings *canvas, Vec3f *canvasVerts, DynamicPaintPainterSettings *paint, Object *canvasOb, Object *paintOb)
{
	DerivedMesh *dm = NULL;
	MVert *mvert;
	MFace *mface;

	if (!paint->dm) {
		printf("DynamicPaint: Invalid paint dm.\n");
		return 0;
	}

	if (paint->flags & MOD_DPAINT_USE_MATERIAL) DynamicPaint_InitMaterialObjects(paintOb);

	{
		BVHTreeFromMesh treeData;
		int yy;
		int tWidth = canvas->surface->w;
		int tHeight = canvas->surface->h;

		int numOfVerts;
		int ii;

		/*
		*	Transform collider vertices to world space
		*	(Faster than transforming per pixel
		*   coordinates and normals to object space)
		*/
		dm = CDDM_copy(paint->dm);
		mvert = dm->getVertArray(dm);
		mface = dm->getFaceArray(dm);
		numOfVerts = dm->getNumVerts(dm);

		for (ii=0; ii<numOfVerts; ii++) {
			mul_m4_v3(paintOb->obmat, mvert[ii].co);
		}

		// Build bvh tree from transformed vertices
		bvhtree_from_mesh_faces(&treeData, dm, 0.0f, 4, 6);


		if(treeData.tree) {


		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < tHeight; yy++)
		{
			int xx;
			for (xx = 0; xx < tWidth; xx++)
			{
				int index = xx+tWidth*yy;
				PaintSurfacePoint *cPoint = (&canvas->surface->point[index]);
				int i = cPoint->index;
				
				// If this canvas point exsist
				if (i >= 0)
				{
					int ss;
					float ssFactor = 0.0f;	// super-sampling factor
					float depth = 0.0f;		// displace depth
					float distRate = -1.0f;	// If it's a proximity hit
											// store distance rate

					float hitCoord[3];		// mid-sample hit coordinate
					int hitFace = -1;		// mid-sample hit face
					short hitQuad;			// mid-sample hit quad status

					// Do supersampling
					for (ss=0; ss<canvas->surface->pixelSamples; ss++) {

						float ray_start[3], ray_dir[3];
						float gaus_factor;
						BVHTreeRayHit hit;
						BVHTreeNearest nearest;
						short hit_found = 0;
						float realPos[3];

						// Supersampling factor
						if (canvas->surface->pixelSamples > 1) {
							gaus_factor = gaussianFactors[ss];
						}
						else {
							gaus_factor = 1.0f;
						}

						// Get current sample position in world coordinates
						interp_v3_v3v3v3(realPos,
										canvasVerts[cPoint->v1].v,
										canvasVerts[cPoint->v2].v,
										canvasVerts[cPoint->v3].v, cPoint->barycentricWeights[ss].v);

						VECCOPY(ray_start, realPos);
						VECCOPY(ray_dir, cPoint->invNorm);


						hit.index = -1;
						hit.dist = 9999;
						nearest.index = -1;
						nearest.dist = paint->paint_distance;

						// Check volume collision
						if (paint->collision == MOD_DPAINT_COL_VOLUME || paint->collision == MOD_DPAINT_COL_VOLDIST)
						if(BLI_bvhtree_ray_cast(treeData.tree, ray_start, ray_dir, 0.0f, &hit, mesh_faces_spherecast_dp, &treeData) != -1)
						{
							// We hit a triangle, now check if collision point normal is facing the point


							//	For optimizion sake, hit point normal isn't calculated in ray cast loop
							int v1=mface[hit.index].v1, v2=mface[hit.index].v2, v3=mface[hit.index].v3, quad=(hit.no[0] == 1.0f);
							float dot;

							if (quad) {v2=mface[hit.index].v3; v3=mface[hit.index].v4;}
							
							// Get hit normal

							normal_tri_v3( hit.no, mvert[v1].co, mvert[v2].co, mvert[v3].co);
							dot = ray_dir[0]*hit.no[0] + ray_dir[1]*hit.no[1] + ray_dir[2]*hit.no[2];

							// If ray and hit normal are facing same direction
							// hit point is inside a closed mesh.
							if (dot>=0)
							{
								// Add factor on supersample filter
								ssFactor += gaus_factor;

								depth += hit.dist * 2.0f;
		
								hit_found = 1;

								/*
								*	Mark first hit info
								*/
								if (hitFace == -1) {
									VECADDFAC(hitCoord, ray_start, ray_dir, hit.dist);
									//copy_v3_v3(hitCoord, hit.co);
									hitQuad = quad;
									hitFace = hit.index;
								}
							}

						}	// end of raycast
						
						// Check proximity collision
						if ((paint->collision == MOD_DPAINT_COL_DIST || paint->collision == MOD_DPAINT_COL_VOLDIST) && (!hit_found))
						{
							float proxDist = -1.0f;
							/*
							*	If pure distance
							*/
							if (!(paint->flags & MOD_DPAINT_PROX_FACEALIGNED)) {
								if (BLI_bvhtree_find_nearest(treeData.tree, ray_start, &nearest, mesh_faces_nearest_point_dp, &treeData) != -1) {
									proxDist = nearest.dist;
								}
							} // non-facealigned
							else {
								negate_v3(ray_dir);
								hit.index = -1;
								hit.dist = paint->paint_distance;
								// Do a face normal directional raycast, and use that distance
								if(BLI_bvhtree_ray_cast(treeData.tree, ray_start, ray_dir, 0.0f, &hit, mesh_faces_spherecast_dp, &treeData) != -1)
								{
									proxDist = hit.dist;
								}
							}

							if (proxDist >= 0.0f) {
								float dist_rate = proxDist / paint->paint_distance;

									// Smooth range or color ramp
									if (paint->proximity_falloff == MOD_DPAINT_PRFALL_SMOOTH ||
										paint->proximity_falloff == MOD_DPAINT_PRFALL_RAMP) {

										if (dist_rate > 1.0f) dist_rate = 1.0f;
										if (dist_rate < 0.0f) dist_rate = 0.0f;
										ssFactor += (1.0f - dist_rate) * gaus_factor;

										if (hitFace == -1) {
											distRate = dist_rate;
										}
									}
									else ssFactor += gaus_factor;

									hit_found = 1;

									if (hitFace == -1) {
										copy_v3_v3(hitCoord, nearest.co);
										hitQuad = (nearest.no[0] == 1.0f);
										hitFace = nearest.index;
									}
							}	// proxDist
						}	// end proximity check



					} // end supersampling


					// if any sample was inside

					if (ssFactor > 0.01f) {

						// apply supersampling results
						if (canvas->surface->pixelSamples > 1) {
							ssFactor /= gaussianTotal;
						}

						cPoint->state = 2;

						if (paint->flags & MOD_DPAINT_DO_PAINT) {

							float paintColor[3];
							float paintAlpha = paint->alpha;
							float paintWetness = paint->wetness * ssFactor;
							float bandres[4];

							paintColor[0] = paint->r;
							paintColor[1] = paint->g;
							paintColor[2] = paint->b;
							
							if (paint->flags & MOD_DPAINT_USE_MATERIAL) DynamicPaint_GetMaterialColor(paintColor, &paintAlpha, paintOb, cPoint->realCoord, hitCoord, hitFace, hitQuad, paint->dm);

							if ((distRate >= 0.0f) && (paint->proximity_falloff == MOD_DPAINT_PRFALL_RAMP) && do_colorband(paint->paint_ramp, distRate, bandres)) {
								if (!(paint->flags & MOD_DPAINT_RAMP_ALPHA)) {
									paintColor[0] = bandres[0];
									paintColor[1] = bandres[1];
									paintColor[2] = bandres[2];
								}
								paintAlpha = bandres[3];
							}

							paintWetness *= paintAlpha;
							paintAlpha = paintAlpha * ssFactor;

							/*
							*	Add paint
							*/
							DynamicPaint_MixPaintColors(cPoint, paint->flags, paintColor, &paintAlpha, &paintWetness);
							
						}

						if (paint->flags & MOD_DPAINT_DO_DISPLACE) {

							float normal_scale, tempNorm[3];
							/*
							*	Calculate normal directional scale of canvas object.
							*	(Displace maps work in object space)
							*/
							MVert *canMvert = NULL;
							canMvert = canvas->dm->getVertArray(canvas->dm);
							normal_tri_v3( tempNorm, canMvert[cPoint->v1].co, canMvert[cPoint->v2].co, canMvert[cPoint->v3].co);
							mul_v3_v3 (tempNorm, canvasOb->size);
							normal_scale = len_v3(tempNorm);
							if (normal_scale<0.01) normal_scale = 0.01;

							depth /= canvas->disp_depth * canvas->surface->pixelSamples * normal_scale;
							if (cPoint->depth < depth) cPoint->depth = depth;

							// Only lift displace if it doesn't have depth
							if (cPoint->depth <= 0.0f && cPoint->depth > depth) cPoint->depth = depth;
						}
					}


				} // end i>0

			} // yy
		} // end of pixel loop
		} // end if tree exists

		// free tree
		free_bvhtree_from_mesh(&treeData);
		dm->release(dm);

	}

	return 1;
}



/*
*	Paint a particle system to canvas
*/
int DynamicPaint_PaintParticles(DynamicPaintCanvasSettings *canvas, ParticleSystem *psys, DynamicPaintPainterSettings *paint, Object *canvasOb)
{
	int yy;
	int tWidth = canvas->surface->w;
	int tHeight = canvas->surface->h;
	ParticleSettings *part=psys->part;
	ParticleData *pa = NULL;

	KDTree *tree;
	int particlesAdded = 0;
	int linvalidParticles = 0;
	int p = 0;


	if (psys->totpart < 1) return 1;

	/*
	*	Build a kd-tree to optimize collision search
	*/
	tree= BLI_kdtree_new(psys->totpart);


	// loop through particles and insert valid ones	to the tree
	for(p=0, pa=psys->particles; p<psys->totpart; p++, pa++)								
				{

					// Proceed only if particle is active
					if(pa->alive == PARS_UNBORN && (part->flag & PART_UNBORN)==0) continue;									
					else if(pa->alive == PARS_DEAD && (part->flag & PART_DIED)==0) continue;									
					else if(pa->flag & PARS_NO_DISP || pa->flag & PARS_UNEXIST) continue;


					// for debug purposes check if any NAN particle proceeds
					// For some reason they get past activity check, this should rule most of them out
					if (isnan(pa->state.co[0]) || isnan(pa->state.co[1]) || isnan(pa->state.co[2])) {linvalidParticles++;continue;}

					BLI_kdtree_insert(tree, p, pa->state.co, NULL);
					particlesAdded++;
					

				}	// particles loop

	if (linvalidParticles) {
		printf("Warning: Invalid particle(s) found!\n");
	}

	if (particlesAdded < 1) {

		BLI_kdtree_free(tree);
		return 1;
	}

	// balance tree
	BLI_kdtree_balance(tree);


	/*
	*	Loop through every pixel
	*/
	#pragma omp parallel for schedule(static)
	for (yy = 0; yy < tHeight; yy++)
	{
		int xx;
		for (xx = 0; xx < tWidth; xx++) {

			int index = xx+tWidth*yy;
			PaintSurfacePoint *cPoint = (&canvas->surface->point[index]);
			int i = cPoint->index;

			// If this canvas point exsist
			if (i >= 0)
			{
					KDTreeNearest nearest;
					int index;

					float distance;
					float solidradius = paint->particle_radius;
					float smooth = paint->particle_smooth;
					float radius = solidradius + smooth;
					float strength;

					// Find nearest particle and get distance to it
					
					index = BLI_kdtree_find_nearest(tree, cPoint->realCoord, NULL, &nearest);
					distance = nearest.dist;

					if (distance > radius) continue;

					// distances inside solid radius should be
					// counted as zero
					distance = (distance-solidradius);
					if (distance<0) distance=0.0f;

					// do smoothness if enabled
					if (smooth) distance/=smooth;

					strength = 1.0f - distance;

					cPoint->state = 2;

					if (paint->flags & MOD_DPAINT_DO_PAINT) {

						float paintAlpha = paint->alpha * strength;
						float paintWetness = paint->wetness * strength;
						float paintColor[3];

						paintColor[0] = paint->r;
						paintColor[1] = paint->g;
						paintColor[2] = paint->b;

						DynamicPaint_MixPaintColors(cPoint, paint->flags, paintColor, &paintAlpha, &paintWetness);

					}

					if (paint->flags & MOD_DPAINT_DO_WETNESS) {
						if (cPoint->wetness < (paint->wetness*strength)) cPoint->wetness = paint->wetness*strength;
					}

					if (paint->flags & MOD_DPAINT_DO_DISPLACE) {
						float sdepth, disp, normal_scale, tempNorm[3];
						/*
						*	Calculate normal directional scale of canvas object.
						*	(Displace maps work in object space)
						*/
						MVert *canMvert = NULL;
						canMvert = canvas->dm->getVertArray(canvas->dm);
						normal_tri_v3( tempNorm, canMvert[cPoint->v1].co, canMvert[cPoint->v2].co, canMvert[cPoint->v3].co);
						mul_v3_v3 (tempNorm, canvasOb->size);
						normal_scale = len_v3(tempNorm);
						if (normal_scale<0.01) normal_scale = 0.01;
						
						sdepth = solidradius - nearest.dist / normal_scale * 2.0f;

						if (sdepth<0.0f) sdepth = 0.0f;
						disp = sdepth / canvas->disp_depth;
						if (cPoint->depth < disp) cPoint->depth = disp;
					}


			}}} // pixel loop and validity check

	BLI_kdtree_free(tree);


	return 1;
}


/*
*	Prepare data required by effects for current frame.
*	Returns number of steps required
*/
static int DynamicPaint_Prepare_EffectStep(DynamicPaintCanvasSettings *canvas, float timescale)
{
	double average_dist;
	float minimum_dist = 9999.0f;
	unsigned int count = 0;
	int steps = 1;
	PaintSurface *surface = canvas->surface;
	int yy, w=surface->w, h=surface->h;

	float fastest_effect, speed_factor;

	/*
	*	Calculate current frame neigh-pixel distances
	*	and average distance between neighbours
	*/
	#pragma omp parallel for schedule(static)
	for (yy = 0; yy < h; yy++)
	{
		int xx;
		for (xx = 0; xx < w; xx++)
		{
			int i;
			PaintSurfacePoint *cPoint = (&surface->point[xx+w*yy]);

			if (cPoint->index != -1)
			for (i=0; i<8; i++) {
				int x,y, index;
				PaintSurfacePoint *tPoint;

				x = xx + neighX[i];
				y = yy + neighY[i];

				index = x+w*y;
				if (cPoint->neighbours[i] != -1) {

					tPoint = (&surface->point[index]);
					cPoint->neighbour_dist[i] = len_v3v3(cPoint->realCoord, tPoint->realCoord);
					if (cPoint->neighbour_dist[i] < minimum_dist) minimum_dist = cPoint->neighbour_dist[i];

					average_dist += cPoint->neighbour_dist[i];
					count++;
				}
			}
		}

	} // end pixel loop
			
	/*
	*	Note: some other method could be better
	*	Right now double precision may cause noticable error on huge
	*	texture sizes (8096*8096*8 = 500 000 000 values)
	*/
	average_dist /= (double)count;

	// Limit minimum distance (used to define substeps)
	// to 1/2 of average
	if (minimum_dist < (average_dist/2.0f)) minimum_dist=average_dist/2.0f;
	//steps = 2;

	fastest_effect = canvas->spread_speed;
	if (canvas->drip_speed > fastest_effect) fastest_effect = canvas->drip_speed;
	if (canvas->shrink_speed > fastest_effect) fastest_effect = canvas->shrink_speed;
	// Number of steps depends on
	// surface res, effect speed and timescale
	speed_factor = (float)surface->w/256.0f * timescale;
	steps = (int)ceil(speed_factor*fastest_effect);

	speed_factor /= (float)steps;



	/*
	*	Now calculate final per pixel speed ratio to neighbour_dist[] array
	*/
	#pragma omp parallel for schedule(static)
	for (yy = 0; yy < h; yy++)
	{
		int xx;
		float dd;
		for (xx = 0; xx < w; xx++)
		{
			int i;
			PaintSurfacePoint *cPoint = (&surface->point[xx+w*yy]);

			if (cPoint->index != -1)
			for (i=0; i<8; i++) {
				if (cPoint->neighbours[i] != -1) {
					dd = (cPoint->neighbour_dist[i]<minimum_dist) ? minimum_dist : cPoint->neighbour_dist[i];
					cPoint->neighbour_dist[i] = minimum_dist / dd * speed_factor;
				}
			}
		}

	} // end pixel loop

	//printf("Average distance is %f, minimum distance is %f\n", average_dist, minimum_dist);

	return steps;
}
/*
*	Clean effect data
*/
void DynamicPaint_Clean_EffectStep(DynamicPaintCanvasSettings *canvas, float timescale)
{
	PaintSurface *surface = canvas->surface;
}
/*
*	Processes active effects.
*/
void DynamicPaint_Do_EffectStep(DynamicPaintCanvasSettings *canvas, PaintSurfacePoint *prevPoint, float timescale)
{
	PaintSurface *surface = canvas->surface;
	int yy, w=surface->w, h=surface->h;

	/*
	*	Spread Effect
	*/
	if (canvas->effect & MOD_DPAINT_EFFECT_DO_SPREAD)  {

		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;
				int i, validPoints = 0;
				float totalAlpha = 0.0f;
				PaintSurfacePoint *cPoint = (&surface->point[index]);	// Current source point
				PaintSurfacePoint *ePoint;							// Effect point to shift values into


				for (i=0; i<8; i++) {
					int nIndex;
					float factor, alphaAdd = 0.0f;

					nIndex = surface->point[index].neighbours[i];

					if (nIndex == -1) continue;

					/*
					*	Find neighbour cells that have higher wetness
					*	and expand it to this cell as well.
					*/
					ePoint = (&prevPoint[nIndex]);
					validPoints++;
					totalAlpha += ePoint->e_alpha;

					if (ePoint->wetness <= cPoint->wetness) continue;
					factor = ePoint->wetness/8 * (ePoint->wetness - cPoint->wetness) * surface->point[index].neighbour_dist[i] * canvas->spread_speed;

					if (ePoint->e_alpha > cPoint->e_alpha) {
						alphaAdd = ePoint->e_alpha/8 * (ePoint->wetness - cPoint->wetness) * surface->point[index].neighbour_dist[i] * canvas->spread_speed;
						//ePoint->e_alpha/8 * (ePoint->e_alpha - cPoint->e_alpha) * surface->point[index].neighbour_dist[i] * canvas->spread_speed;
					}


					if (cPoint->e_alpha) {
						float invFactor = 1.0f - factor;
						cPoint->e_color[0] = cPoint->e_color[0]*invFactor + ePoint->e_color[0]*factor;
						cPoint->e_color[1] = cPoint->e_color[1]*invFactor + ePoint->e_color[1]*factor;
						cPoint->e_color[2] = cPoint->e_color[2]*invFactor + ePoint->e_color[2]*factor;
							
						cPoint->e_alpha += alphaAdd;
						cPoint->wetness += factor;
					}
					else {
						// If there is no existing paint, just replace the color
						cPoint->e_color[0] = ePoint->e_color[0];
						cPoint->e_color[1] = ePoint->e_color[1];
						cPoint->e_color[2] = ePoint->e_color[2];

						cPoint->e_alpha += alphaAdd;
						cPoint->wetness += factor;
					}

				}

				// For antialiasing sake, don't let alpha go higher than average alpha of neighbours
				if (validPoints && (cPoint->e_alpha > (totalAlpha/validPoints+0.25f))) {
					cPoint->e_alpha = (totalAlpha/validPoints+0.25f);
					if (cPoint->e_alpha>1.0f) cPoint->e_alpha = 1.0f;
				}

			}
		} // end pixel loop
	} // end spread effect






	/*
	*	Spread Effect
	*/
	if (canvas->effect & MOD_DPAINT_EFFECT_DO_DRIP) 
	{	// Dripping

		/*
		*	Test for generating randomized dripping effect:
		*	Generates random areas of higher wetness to
		*	the surface and thus makes some parts drip faster.
		*
		*	Not properly implemented yet.
		*/
#if 0
		{
				KDTree *tree;
				int id = 0;
				int max_count = 4096;
				int random_factor =

				/*
				*	Build a kd-tree to optimize search
				*/
				tree= BLI_kdtree_new(4096);	// just use
				srand ( 32 );

				/*
				*	Loop through every pixel and insert drip points
				*	to random point real world coordinates
				*	(Non threadable?)
				*/
				/*for (yy = 0; yy < surface->h; yy++)
				{
					int xx;
					for (xx = 0; xx < surface->w; xx++) {

						int index = xx+surface->w*yy;
						PaintSurfacePoint *cPoint = (&canvas->surface->point[index]);
						int i = cPoint->index;

						// If this canvas point exsist
						if (i >= 0)
						{
								if ((rand() % 300) == 0) {

									// Add a flow point on every 1.000th pixel

									BLI_kdtree_insert(tree, id, cPoint->realCoord, NULL);
									id++;
								}

								//index = BLI_kdtree_find_nearest(tree, cPoint->realCoord, NULL, &nearest);
								//distance = nearest.dist;


				}}} // pixel loop and validity check

					// balance tree
				BLI_kdtree_balance(tree);*/

				/*
				*	Loop through every pixel and double wetness on new paint areas
				*	near defined drip points
				*/
				/*#pragma omp parallel for schedule(static)
				for (yy = 0; yy < surface->h; yy++)
				{
					int xx;
					for (xx = 0; xx < surface->w; xx++) {

						int index = xx+surface->w*yy;
						PaintSurfacePoint *cPoint = (&canvas->surface->point[index]);
						int i = cPoint->index;

						// If this canvas point exsist
						if (i >= 0 && cPoint->state == 2)
						{
								KDTreeNearest nearest;
								int index;

								float distance;

								index = BLI_kdtree_find_nearest(tree, cPoint->realCoord, NULL, &nearest);
								distance = nearest.dist;

								if (distance < 0.20) {
									float strength = (1.0f - (distance / 0.20));
									if (strength > 0.8f) strength = 0.8f;
									strength *= 1.6f;
									strength += 1.0f;

									cPoint->wetness *= strength;
									if (cPoint->wetness > 2.0) cPoint->wetness = 2.0f;
									cPoint->e_color[2] = strength - 1.0f;
								}


				}}} // pixel loop and validity check

				BLI_kdtree_free(tree);*/
		}	// end dripping random gene (alpha)
#endif




		memcpy(prevPoint, surface->point, surface->w*surface->h*sizeof(struct PaintSurfacePoint));

		//#pragma omp parallel for schedule(static)
		for (yy = 1; yy < h-1; yy++)
		{
			int xx;
			for (xx = 1; xx < w-1; xx++)
			{
				int index = xx+w*yy;
				int validPoints = 0;
				float totalFactor = 0.0f, totalAlpha = 0.0f;
				float factor = 0.0f, drip_strength = 0.0f;
				PaintSurfacePoint *cPoint = (&prevPoint[index]);	// Current source point
				PaintSurfacePoint *dPoint = (&surface->point[index]);	// Current dest point
				PaintSurfacePoint *ePoint;							// Effect point to shift values into
				PaintSurfacePoint *ePoint2;							// Effect point to shift values into


				int nIndex;

				// Get dir modulo PI/2, and adjust direction to match halfway of two dir pixels by substracting 3/4PI

				// Because the neighbour pixels go in 45 degree slices
				// Get drip factor in that range
				float dirMod = fmod(cPoint->gravity_dir,0.785398163);
				float stFac=dirMod/0.785398163;
				float facTotal;

				// Get neighpixel index
				int neighPixel = (int)floor(cPoint->gravity_dir/6.2831853f * 8.0f);
				int neighPixel2;

				if (neighPixel > 7) neighPixel = 7;	// Shouldn't happen but just in case
				if (neighPixel < 0) neighPixel = 0;


				// Make sure neighbour exists
				nIndex = surface->point[index].neighbours[neighPixel];
				if (nIndex == -1) continue;
				ePoint = (&surface->point[nIndex]);

				// Do some adjustments to dripping speed
				drip_strength = cPoint->wetness - 0.2f;
				if (drip_strength < 0) continue;



				/*
				*	Second neighbour point
				*/
				neighPixel2 = neighPixel + 1;
				if (neighPixel2 > 7) neighPixel2 = 0;
				nIndex = surface->point[index].neighbours[neighPixel2];
				if (nIndex == -1) continue;
				ePoint2 = (&surface->point[nIndex]);

				factor = 0.3 * canvas->drip_speed;
				if (drip_strength < 2.5f) factor *= drip_strength;

				// Add values to pixel below
				ePoint->e_color[0] = cPoint->e_color[0];
				ePoint->e_color[1] = cPoint->e_color[1];
				ePoint->e_color[2] = cPoint->e_color[2];

				ePoint->e_alpha += cPoint->e_alpha * factor * surface->point[index].neighbour_dist[neighPixel] * (1.0f - stFac);
				ePoint->wetness += cPoint->wetness * factor * surface->point[index].neighbour_dist[neighPixel] * (1.0f - stFac);

				if (stFac > 0.0f) {
					ePoint2->e_color[0] = cPoint->e_color[0];
					ePoint2->e_color[1] = cPoint->e_color[1];
					ePoint2->e_color[2] = cPoint->e_color[2];
				}

				ePoint2->e_alpha += cPoint->e_alpha * factor * surface->point[index].neighbour_dist[neighPixel2] * stFac;
				ePoint2->wetness += cPoint->wetness * factor * surface->point[index].neighbour_dist[neighPixel2] * stFac;

				// Decrease values on this one
				//dPoint->e_color[0] -= cPoint->e_color[0] * factor;
				//dPoint->e_color[1] -= cPoint->e_color[1] * factor;
				//dPoint->e_color[2] -= cPoint->e_color[2] * factor;

				//dPoint->e_alpha -= cPoint->e_alpha * factor;
				facTotal = surface->point[index].neighbour_dist[neighPixel] * (1.0f - stFac) + surface->point[index].neighbour_dist[neighPixel2] * stFac;
				dPoint->wetness -= cPoint->wetness * factor * facTotal;

			}
		} // end pixel loop


		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;

				PaintSurfacePoint *cPoint = (&surface->point[index]);	// Current source point

				if (cPoint->e_alpha > 1.0f) cPoint->e_alpha=1.0f;
				if (cPoint->wetness > 3.5f) cPoint->wetness=3.5f;

				if (cPoint->e_alpha < 0.0f) cPoint->e_alpha=0.0f;
				if (cPoint->wetness < 0.0f) cPoint->wetness=0.0f;

			}
		} // end pixel loop
	} // end dripping effect



	/*
	*	Shrink Effect
	*/
	if (canvas->effect & MOD_DPAINT_EFFECT_DO_SHRINK)  {

		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < h; yy++)
		{
			int xx;
			for (xx = 0; xx < w; xx++)
			{
				int index = xx+w*yy;
				int i, validPoints = 0;
				float totalAlpha = 0.0f;
				PaintSurfacePoint *cPoint = (&surface->point[index]);	// Current source point
				PaintSurfacePoint *ePoint;							// Effect point to shift values into


				for (i=0; i<8; i++) {
					int nIndex;
					float factor, e_factor, w_factor;

					nIndex = surface->point[index].neighbours[i];

					if (nIndex == -1) continue;

					/*
					*	Find neighbour cells that have lower alpha
					*	and decrease this point alpha towards that level.
					*/
					ePoint = (&prevPoint[nIndex]);
					validPoints++;
					totalAlpha += ePoint->e_alpha;

					if (cPoint->alpha <= 0.0f && cPoint->e_alpha <= 0.0f && cPoint->wetness <= 0.0f) continue;
					factor = (1.0f - ePoint->alpha)/8 * (cPoint->alpha - ePoint->alpha) * surface->point[index].neighbour_dist[i] * canvas->shrink_speed;
					if (factor < 0.0f) factor = 0.0f;

					e_factor = (1.0f - ePoint->e_alpha)/8 * (cPoint->e_alpha - ePoint->e_alpha) * surface->point[index].neighbour_dist[i] * canvas->shrink_speed;
					if (e_factor < 0.0f) e_factor = 0.0f;

					w_factor = (1.0f - ePoint->wetness)/8 * (cPoint->wetness - ePoint->wetness) * surface->point[index].neighbour_dist[i] * canvas->shrink_speed;
					if (w_factor < 0.0f) w_factor = 0.0f;


					if (factor) {
						cPoint->alpha -= factor;
						if (cPoint->alpha < 0.0f) cPoint->alpha = 0.0f;
						cPoint->wetness -= factor;

					}
					else {
						cPoint->e_alpha -= e_factor;
						if (cPoint->e_alpha < 0.0f) cPoint->e_alpha = 0.0f;
						cPoint->wetness -= w_factor;
						if (cPoint->wetness < 0.0f) cPoint->wetness = 0.0f;
					}

				}

			}
		} // end pixel loop
	} // end spread effect
}



/*
*	Do Paint Step. Paints scene objects of current frame.
*/
int DynamicPaint_DoStep(Scene *scene, Object *ob, DynamicPaintModifierData *pmd, float timescale)
{
	DynamicPaintCanvasSettings *canvas = pmd->canvas;
	PaintSurface *surface = canvas->surface;
	Base *base = NULL;
	MVert *mvert = NULL;
	MFace *mface = NULL;
	MTFace *tface = NULL;
	DerivedMesh *dm = canvas->dm;

	int w, h;
	int yy;

	int canvasNumOfVerts;
	struct Vec3f *canvasVerts;
	int canvasNumOfFaces;
	struct FaceAdv *canvasInvNormals;

	if (!dm) {
		printf("DynamicPaint : Invalid canvas mesh\n");
		return 0;
	}

	mvert = dm->getVertArray(dm);
	mface = dm->getFaceArray(dm);
	tface = DM_get_face_data_layer(dm, CD_MTFACE);
	canvasNumOfFaces = dm->getNumFaces(dm);

	w = canvas->surface->w;
	h = canvas->surface->h;

	/*
	*	Make a transformed copy of canvas derived mesh vertices to avoid recalculation.
	*/

	canvasNumOfVerts = dm->getNumVerts(dm);
	canvasVerts =  (struct Vec3f *) MEM_mapallocN(canvasNumOfVerts*sizeof(struct Vec3f), "Dynamic Paint Transformed canvas");

	#pragma omp parallel for schedule(static)
	for (yy=0; yy<canvasNumOfVerts; yy++) {
		VECCOPY(canvasVerts[yy].v, mvert[yy].co);
		// Multiply coordinates by canvas matrix
		mul_m4_v3(ob->obmat, canvasVerts[yy].v);
	}


	/*
	*	Calculate temp per face normals using transformed coordinates.
	*	(To not have to calculate same normal for millions of pixels)
	*/
	canvasInvNormals =  (struct FaceAdv *) MEM_mapallocN(canvasNumOfFaces*sizeof(struct FaceAdv), "Dynamic Paint canvas normals");

	#pragma omp parallel for schedule(static)
	for (yy=0; yy<canvasNumOfFaces; yy++) {

		if (mface[yy].flag & ME_SMOOTH) continue;	// Only calculate flat faces
		/*
		*	Transformed normal
		*/
		normal_tri_v3( canvasInvNormals[yy].no, canvasVerts[mface[yy].v3].v, canvasVerts[mface[yy].v2].v, canvasVerts[mface[yy].v1].v);
		if (mface[yy].v4) normal_tri_v3( canvasInvNormals[yy].no_q, canvasVerts[mface[yy].v4].v, canvasVerts[mface[yy].v3].v, canvasVerts[mface[yy].v1].v);
	}


	/*
	*	Do per pixel stuff like dissolve effect and calculate pixel center coordinates
	*/
	#pragma omp parallel for schedule(static)
	for (yy = 0; yy < h; yy++)
	{
		int xx;
		float n1[3], n2[3], n3[3];
		for (xx = 0; xx < w; xx++)
		{
			int index = xx+w*yy;
			short i, quad;


			i = canvas->surface->point[index].index;
			quad = canvas->surface->point[index].quad;


			// if current pixel has 3D coordinates
			if (i >= 0)
			{
				PaintSurfacePoint *cPoint = (&canvas->surface->point[index]);


				/*
				*  Do dissolve / drying
				*/

				if (cPoint->wetness > 0.0f) {

					// Every dryinh step Blends wet paint to the background.
					float invAlpha = 1.0f - cPoint->e_alpha;
					cPoint->color[0] = cPoint->color[0]*invAlpha + cPoint->e_color[0]*cPoint->e_alpha;
					cPoint->color[1] = cPoint->color[1]*invAlpha + cPoint->e_color[1]*cPoint->e_alpha;
					cPoint->color[2] = cPoint->color[2]*invAlpha + cPoint->e_color[2]*cPoint->e_alpha;

					cPoint->state = 1;

					// only increase alpha if wet paint has higher
					if (cPoint->e_alpha > cPoint->alpha) cPoint->alpha = cPoint->e_alpha;

					// Now dry it ;o
					if (canvas->flags & MOD_DPAINT_DRY_LOG) cPoint->wetness *= 1.0f - (1.0 / (canvas->dry_speed/timescale));
					else cPoint->wetness -= 1.0f/canvas->dry_speed*timescale;
				}
				if (cPoint->wetness <= 0.0f) {
					cPoint->wetness = 0.0f;
					cPoint->e_alpha = 0.0f;
					cPoint->state = 0;
				}

				if (canvas->flags & (MOD_DPAINT_DISSOLVE)) {

					cPoint->alpha -= 1.0f/canvas->diss_speed*timescale;
					if (cPoint->alpha < 0.0f) cPoint->alpha = 0.0f;

					cPoint->e_alpha -= 1.0f/canvas->diss_speed*timescale;
					if (cPoint->e_alpha < 0.0f) cPoint->e_alpha = 0.0f;
				}

				if (canvas->flags & (MOD_DPAINT_FLATTEN)) {

					cPoint->depth -= 1.0f/canvas->dflat_speed*timescale;
					if (cPoint->depth < 0.0f) cPoint->depth = 0.0f;
				}





				/*
				*	Calculate current 3D-position of each texture pixel
				*/

				// Save pixel real coord to canvas data
				interp_v3_v3v3v3(	cPoint->realCoord,
					canvasVerts[cPoint->v1].v,
					canvasVerts[cPoint->v2].v,
					canvasVerts[cPoint->v3].v, cPoint->barycentricWeights[0].v);

				// Calculate pixel normal
				if(mface[cPoint->index].flag & ME_SMOOTH) {
					normal_short_to_float_v3(n1, mvert[cPoint->v1].no);
					normal_short_to_float_v3(n2, mvert[cPoint->v2].no);
					normal_short_to_float_v3(n3, mvert[cPoint->v3].no);

					interp_v3_v3v3v3(	cPoint->invNorm,
						n1, n2, n3, cPoint->barycentricWeights[0].v);
					normalize_v3(cPoint->invNorm);
				}
				else {
					if (cPoint->quad) {VECCOPY(cPoint->invNorm, canvasInvNormals[cPoint->index].no_q);}
					else {VECCOPY(cPoint->invNorm, canvasInvNormals[cPoint->index].no);}
				}


				if (canvas->effect) {
					/*
					*	Get current gravity direction of pixel in UV space.
					*/
					cPoint->gravity_dir = 0.0f;
					cPoint->gravity_rate = 1.0f - abs(cPoint->invNorm[2]);

					// World gravity direction is negative z-axis
					if (cPoint->gravity_rate > 0.01)
					{
						float uv1[3], uv2[3], uv3[3], unormal[3];
						int v2i, v3i;

						v2i = (cPoint->quad) ? 2 : 1;
						v3i = (cPoint->quad) ? 3 : 2;

						/*
						*	Apply z-coodrinate (gravity dir) to uv-vertices
						*/
						uv1[0] = tface[cPoint->index].uv[0][0];
						uv1[1] = tface[cPoint->index].uv[0][1];
						uv1[2] = canvasVerts[cPoint->v1].v[2];

						uv2[0] = tface[cPoint->index].uv[v2i][0];
						uv2[1] = tface[cPoint->index].uv[v2i][1];
						uv2[2] = canvasVerts[cPoint->v2].v[2];

						uv3[0] = tface[cPoint->index].uv[v3i][0];
						uv3[1] = tface[cPoint->index].uv[v3i][1];
						uv3[2] = canvasVerts[cPoint->v3].v[2];

						/*
						*	Calculate a new normal vector for that generated face
						*/
						normal_tri_v3( unormal, uv1, uv2, uv3);

						// Normalize u and v part of it
						normalize_v2(unormal);

						/*
						*	And use direction of that 2d vector as gravity direction in uv space
						*	(Convert from -pi->pi to 0->2pi)
						*/
						cPoint->gravity_dir = fmod((atan2(unormal[1], unormal[0])+6.28318531f), 6.28318531f);

					}
				} // end if (effect)

			} // end i>0


		}	// yy
	} // end of loop xx

	MEM_freeN(canvasInvNormals);

	/*
	*	Do collision checks for every painter
	*/
	{
		Object *otherobj = NULL;
		ModifierData *md = NULL;

		base = scene->base.first;

		while(base)
		{

			otherobj = base->object;

			if(!otherobj)					
			{												
				base= base->next;						
				continue;			
			}
			base= base->next;

			md = modifiers_findByType(otherobj, eModifierType_DynamicPaint);

			// check if target has an active mp modifier
			if(md && md->mode & (eModifierMode_Realtime | eModifierMode_Render))					
			{
				DynamicPaintModifierData *pmd2 = (DynamicPaintModifierData *)md;


				// Make sure we're dealing with a painter
				if((pmd2->type & MOD_DYNAMICPAINT_TYPE_PAINT) && pmd2->paint)
				{

					DynamicPaintPainterSettings *paint = pmd2->paint;

					/*
					*	Make sure at least one influence is enabled
					*/
					if (!(paint->flags & MOD_DPAINT_DO_PAINT || paint->flags & MOD_DPAINT_DO_DISPLACE)) continue;


					// Check if painter has a particle system selected
					if (paint->collision == MOD_DPAINT_COL_PSYS)
					{

						if (paint && paint->psys && paint->psys->part && paint->psys->part->type==PART_EMITTER)
						if (psys_check_enabled(otherobj, paint->psys)) {
							/*
							*	Paint the particle system
							*/
							DynamicPaint_PaintParticles(canvas, paint->psys, paint, ob);
						}
										
					}							
					else {

						/*
						*	Paint the mesh
						*/
						DynamicPaint_PaintMesh(canvas, canvasVerts, paint, ob, otherobj);
					}


					} // end of collision check (Is valid coll modifier)

			}
			}
		}

		// Free per frame canvas data
		MEM_freeN(canvasVerts);








		/*
		*	EFFECT HANDLING
		*/

		if (canvas->effect)
		{
			unsigned int steps = 1, s;

			// Allocate memory for surface previous points
			PaintSurfacePoint *prevPoint = (struct PaintSurfacePoint *) MEM_mapallocN(surface->w*surface->h*sizeof(struct PaintSurfacePoint), "PaintSurfaceDataTemp");


			steps = DynamicPaint_Prepare_EffectStep(canvas, timescale);


			/*
			*	Do Effects steps
			*/
			for (s = 0; s < steps; s++)
			{
				// Copy current surface to previous array
				memcpy(prevPoint, surface->point, surface->w*surface->h*sizeof(struct PaintSurfacePoint));

				DynamicPaint_Do_EffectStep(canvas, prevPoint, timescale);

			}

			MEM_freeN(prevPoint);

			DynamicPaint_Clean_EffectStep(canvas, timescale);

		} // end effects

		return 1;
}


/*
*	Outputs an image of canvas surface
*/
void dynamic_paint_output_image(struct DynamicPaintCanvasSettings *canvas, char* filename, short type, short source)
{
		int yy;
		ImBuf* mhImgB = NULL;

		int tWidth, tHeight;
		PaintSurface *surface = canvas->surface;


		if (!surface) {printf("Canvas save failed: Invalid surface.\n");return;}

		tWidth = surface->w;
		tHeight = surface->h;

		// Init image buffer
		mhImgB = IMB_allocImBuf(tWidth, tHeight, 32, IB_rectfloat, 0);

		#pragma omp parallel for schedule(static)
		for (yy = 0; yy < tHeight; yy++)
		{
			int xx;
			for (xx = 0; xx < tWidth; xx++)
			{
				int pos=(xx+tWidth*yy)*4;	// image buffer position
				int index=(xx+tWidth*yy);	// surface point
				PaintSurfacePoint *cPoint = (&surface->point[index]);

				// if point has no index but has a neighbour, use it
				//if (surface->point[index].neighbour_index != -1) index = surface->point[index].neighbour_index;
				

				// Set values of preferred type
				if (source == MPOUTPUT_WET) {
					float value = (cPoint->wetness > 1.0f) ? 1.0f : cPoint->wetness;
					mhImgB->rect_float[pos]=value;
					mhImgB->rect_float[pos+1]=value;
					mhImgB->rect_float[pos+2]=value;
					mhImgB->rect_float[pos+3]=1.0f;
				}
				else if (source == MPOUTPUT_DISPLACE) {

					float depth = cPoint->depth;

					if (canvas->disp_type & MOD_DPAINT_DISP_DISPLACE) {
						depth = (1.0f - depth) / 2;
					}

					mhImgB->rect_float[pos]=depth;
					mhImgB->rect_float[pos+1]=depth;
					mhImgB->rect_float[pos+2]=depth;
					mhImgB->rect_float[pos+3]=1.0f;
				}else {	// otherwise output MPOUTPUT_PAINT

					float invAlpha = 1.0f - cPoint->e_alpha;

					// If base layer already has a color, blend it
					if (cPoint->alpha) {
						mhImgB->rect_float[pos]   = cPoint->color[0] * invAlpha + cPoint->e_color[0] * cPoint->e_alpha;
						mhImgB->rect_float[pos+1] = cPoint->color[1] * invAlpha + cPoint->e_color[1] * cPoint->e_alpha;
						mhImgB->rect_float[pos+2] = cPoint->color[2] * invAlpha + cPoint->e_color[2] * cPoint->e_alpha;
					}
					else {
						// Else use effect layer color
						mhImgB->rect_float[pos]   = cPoint->e_color[0];
						mhImgB->rect_float[pos+1] = cPoint->e_color[1];
						mhImgB->rect_float[pos+2] = cPoint->e_color[2];
					}

					// Set alpha to the highest
					mhImgB->rect_float[pos+3] = (cPoint->e_alpha > cPoint->alpha) ? cPoint->e_alpha : cPoint->alpha;

					// Multiply color by alpha if set
					if (canvas->flags & MOD_DPAINT_MULALPHA) {
						mhImgB->rect_float[pos]   *= mhImgB->rect_float[pos+3];
						mhImgB->rect_float[pos+1] *= mhImgB->rect_float[pos+3];
						mhImgB->rect_float[pos+2] *= mhImgB->rect_float[pos+3];
					}

				}

			}
		}

		// Save image buffer
		if (type == MPOUTPUT_JPEG) {	// JPEG
			mhImgB->ftype= JPG|95;
			IMB_rect_from_float(mhImgB);
			imb_savejpeg(mhImgB, filename, IB_rectfloat); }
#ifdef WITH_OPENEXR
		else if (type == MPOUTPUT_OPENEXR) {	// OpenEXR 16-bit half float
			mhImgB->ftype = OPENEXR_HALF | OPENEXR_COMPRESS;
			IMB_rect_from_float(mhImgB);
			imb_save_openexr(mhImgB, filename, IB_rectfloat); }
#endif
		else {	// MPOUTPUT_PNG
			mhImgB->ftype= PNG|95;
			IMB_rect_from_float(mhImgB);
			imb_savepng(mhImgB, filename, IB_rectfloat); }


		IMB_freeImBuf(mhImgB);
}


/*
*	Do actual bake operation.
*	returns 0 on failture
*/
int DynamicPaint_Bake(bContext *C, struct DynamicPaintModifierData *pmd, Object *cObject)
{

	DynamicPaintCanvasSettings *canvas;
	Scene *scene= CTX_data_scene(C);
	wmWindow *win = CTX_wm_window(C);
	int frame = 1;
	int frames;

	canvas = pmd->canvas;
	if (!canvas) {printf("DynamicPaint bake failed: Invalid canvas.\n"); return 0;}

	frames = canvas->end_frame - canvas->start_frame + 1;
	if (frames <= 0) {printf("DynamicPaint: No frames to bake.\n"); return 0;}

	/*
	*	Set frame to start point (also inits modifier data) and create surface
	*/
	frame = canvas->start_frame;
	scene->r.cfra = (int)frame;
	ED_update_for_newframe(C, 1);

	// Init surface
	if (!dynamicPaint_createCanvasSurface(canvas)) return 0;

	/*
	*	Loop through scene frames
	*/
	for (frame=canvas->start_frame; frame<=canvas->end_frame; frame++)
	{
		float timescale = 1.0f / (canvas->substeps+1);
		int st;
		float progress = (frame - 1) / (float)frames * 100;

		// If user requested stop (esc), quit baking
		if (blender_test_break()) return 0;

		// Update progress bar cursor
		WM_timecursor(win, (int)progress);

		printf("DynamicPaint: Baking frame %i\n", frame);


		/*
		*	Do calculations for every substep
		*	Note: these have to be from previous frame
		*/
		if (frame != canvas->start_frame) {
			for (st = 1; st <= canvas->substeps; st++)
			{
				float subframe = ((float) st) / (canvas->substeps+1);
				// Update frame if we have proceed

				scene->r.cfra = (int)frame - 1;
				scene->r.subframe = subframe;
				ED_update_for_newframe(C, 1);

				if (!DynamicPaint_DoStep(scene, cObject, pmd, timescale)) return 0;
			}

			/*
			*	Calculate actual frame data
			*/
			scene->r.cfra = (int)frame;
			scene->r.subframe = 0.0f;
			ED_update_for_newframe(C, 1);

		}

		if (!DynamicPaint_DoStep(scene, cObject, pmd, timescale)) return 0;



		/*
		*	Just in case, check if any output is enabled
		*	Don't cancel the bake as user may have keyframed outputs
		*/
		if (!(canvas->output & MOD_DPAINT_OUT_PAINT || canvas->output & MOD_DPAINT_OUT_WET || canvas->output & MOD_DPAINT_OUT_DISP)) {
			printf("Skipping output.\n", frame);
			continue;
		}
		/*
		*	Save output
		*/
		{
		char filename[250];
		char pad[4];
		char wet[4];
		char disp[4];

		// Add frame number padding
		if (frame<10) sprintf(pad,"000");
		else if (frame<100) sprintf(pad,"00");
		else if (frame<1000) sprintf(pad,"0");
		else pad[0] = '\0';

		// Check if paint and wet map filename is same and fix if necessary
		if (!strcmp(canvas->paint_output_path, canvas->wet_output_path)) sprintf(wet,"wet");
		else wet[0] = '\0';
		// same for displacement map
		if (!strcmp(canvas->paint_output_path, canvas->displace_output_path)) sprintf(disp,"disp");
		else if (!strcmp(canvas->wet_output_path, canvas->displace_output_path)) sprintf(disp,"disp");
		else disp[0] = '\0';


		// color map
		if (canvas->output & MOD_DPAINT_OUT_PAINT) {
			sprintf(filename, "%s%s%i.png", canvas->paint_output_path, pad, (int)frame);
			BLI_path_abs(filename, G.sce);
			BLI_make_existing_file(filename);

			printf("Saving color map : %s\n", filename);
			dynamic_paint_output_image(canvas, filename, MPOUTPUT_PNG, MPOUTPUT_PAINT);
		}

		// wetness map
		if (canvas->output & MOD_DPAINT_OUT_WET) {
			sprintf(filename, "%s%s%s%i.png", canvas->wet_output_path, wet, pad, (int)frame);
			BLI_path_abs(filename, G.sce);
			BLI_make_existing_file(filename);

			printf("Saving wet map to : %s\n", filename);
			dynamic_paint_output_image(canvas, filename, MPOUTPUT_PNG, MPOUTPUT_WET);
		}

		// displacement map
		if (canvas->output & MOD_DPAINT_OUT_DISP) {
			// OpenEXR or PNG
			int format = (canvas->disp_format & MOD_DPAINT_DISPFOR_OPENEXR) ? MPOUTPUT_OPENEXR : MPOUTPUT_PNG;
			char ext[4];
			if (canvas->disp_format & MOD_DPAINT_DISPFOR_OPENEXR) sprintf(ext,"exr"); else sprintf(ext,"png");

			sprintf(filename, "%s%s%s%i.%s", canvas->displace_output_path, disp, pad, (int)frame, ext);
			BLI_path_abs(filename, G.sce);
			BLI_make_existing_file(filename);

			printf("Saving displacement map to : %s\n", filename);
			dynamic_paint_output_image(canvas, filename, format, MPOUTPUT_DISPLACE);
		}

		}
	}

	return 1;
}


/*
*	Updates baking status for every paint object
*	returns 0 if no paint objects found.
*/
int DynamicPaint_PainterSetBaking(Scene *scene, short baking) {
	Object *otherobj = NULL;
	ModifierData *md = NULL;
	int count = 0;

	Base *base = scene->base.first;

	while(base)
	{

		otherobj = base->object;

		if(!otherobj)					
		{												
			base= base->next;						
			continue;			
		}

		md = modifiers_findByType(otherobj, eModifierType_DynamicPaint);

		// check if target has an active mp modifier
		if(md && md->mode & (eModifierMode_Realtime | eModifierMode_Render))					
		{
			DynamicPaintModifierData *pmd2 = (DynamicPaintModifierData *)md;
			// Make sure we're dealing with a painter
			if((pmd2->type & MOD_DYNAMICPAINT_TYPE_PAINT) && pmd2->paint)
			{
				pmd2->baking = baking;
				count++;
			}
		}

		base= base->next;
	}

	if (count) return 1;
	
	return 0;
}


/*
*	An operator call to bake dynamic paint simulation for active object
*/
int dynamic_paint_bake_all(bContext *C, Object *cObject) {

	DynamicPaintModifierData *pmd = NULL;
	Scene *scene= CTX_data_scene(C);
	int status = 0;
	double timer = PIL_check_seconds_timer();
	char timestr[30];

	/*
	*	Get modifier data
	*/
	pmd = (DynamicPaintModifierData *)modifiers_findByType(cObject, eModifierType_DynamicPaint);
	if (!pmd) {printf("DynamicPaint bake failed: No Dynamic Paint modifier found.\n"); return 0;}


	/*
	*	Set state to baking and init surface
	*/
	pmd->baking = 1;
	if (!DynamicPaint_PainterSetBaking(scene, 1)) {printf("DynamicPaint bake failed: No paint objects.\n"); return 0;}
	G.afbreek= 0;	/* reset blender_test_break*/

	status = DynamicPaint_Bake(C, pmd, cObject);

	/*
	*	Clean bake stuff
	*/
	pmd->baking = 0;
	DynamicPaint_PainterSetBaking(scene, 0);
	dynamicPaint_cleanCanvasSurface(pmd->canvas);
	// Restore cursor back to normal
	WM_cursor_restore(CTX_wm_window(C));

	/*
	*	Report for ended bake
	*/
	{
		// Format time string
		double time = PIL_check_seconds_timer() - timer;
		int tmp_val;
		timestr[0] = '\0';

		// days (just in case someone actually has a very slow pc)
		tmp_val = (int)floor(time / 86400.0f);
		if (tmp_val > 0) sprintf(timestr, "%i Days - ", tmp_val);
		// hours
		time -= 86400.0f * tmp_val;
		tmp_val = (int)floor(time / 3600.0f);
		if (tmp_val > 0) sprintf(timestr, "%s%i h ", timestr, tmp_val);
		// minutes
		time -= 3600.0f * tmp_val;
		tmp_val = (int)floor(time / 60.0f);
		if (tmp_val > 0) sprintf(timestr, "%s%i min ", timestr, tmp_val);
		// seconds
		time -= 60.0f * tmp_val;
		tmp_val = (int)ceil(time);
		sprintf(timestr, "%s%i s", timestr, tmp_val);
	}
	if (status) printf("Baking Complete! (Time: %s)\n", timestr);
	else  printf("Baking Cancelled!\n");

	return status;
}


/***************************** Operators ******************************/

static int dynamicpaint_bake_exec(bContext *C, wmOperator *op)
{
	Object *ob= CTX_data_active_object(C);

	if(!dynamic_paint_bake_all(C, ob))
		return OPERATOR_CANCELLED;

	return OPERATOR_FINISHED;
}

void DPAINT_OT_bake(wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Dynamic Paint Bake";
	ot->description= "Bake dynamic paint";
	ot->idname= "DPAINT_OT_bake";
	
	/* api callbacks */
	ot->exec= dynamicpaint_bake_exec;
	ot->poll= ED_operator_object_active_editable;
}


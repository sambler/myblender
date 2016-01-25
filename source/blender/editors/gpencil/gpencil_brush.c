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
 * The Original Code is Copyright (C) 2015, Blender Foundation
 * This is a new part of Blender
 *
 * Contributor(s): Joshua Leung
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 * Brush based operators for editing Grease Pencil strokes
 */

/** \file blender/editors/gpencil/gpencil_brush.c
 *  \ingroup edgpencil
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_ghash.h"
#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_utildefines.h"

#include "BLT_translation.h"

#include "DNA_scene_types.h"
#include "DNA_screen_types.h"
#include "DNA_space_types.h"
#include "DNA_view3d_types.h"
#include "DNA_gpencil_types.h"

#include "BKE_context.h"
#include "BKE_global.h"
#include "BKE_gpencil.h"
#include "BKE_library.h"
#include "BKE_report.h"
#include "BKE_screen.h"

#include "UI_interface.h"

#include "WM_api.h"
#include "WM_types.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "UI_view2d.h"

#include "ED_gpencil.h"
#include "ED_screen.h"
#include "ED_view3d.h"

#include "BIF_gl.h"
#include "BIF_glutil.h"

#include "gpencil_intern.h"

/* ************************************************ */
/* General Brush Editing Context */

/* Context for brush operators */
typedef struct tGP_BrushEditData {
	/* Current editor/region/etc. */
	/* NOTE: This stuff is mainly needed to handle 3D view projection stuff... */
	Scene *scene;
	
	ScrArea *sa;
	ARegion *ar;
	
	/* Current GPencil datablock */
	bGPdata *gpd;
	
	/* Brush Settings */
	GP_BrushEdit_Settings *settings;
	GP_EditBrush_Data *brush;
	
	eGP_EditBrush_Types brush_type;
	eGP_EditBrush_Flag  flag;
	
	/* Space Conversion Data */
	GP_SpaceConversion gsc;
	
	
	/* Is the brush currently painting? */
	bool is_painting;
	
	/* Start of new sculpt stroke */
	bool first;
	
	/* Current frame */
	int cfra;
	
	
	/* Brush Runtime Data: */
	/* - position and pressure
	 * - the *_prev variants are the previous values
	 */
	int   mval[2], mval_prev[2];
	float pressure, pressure_prev;
	
	/* - effect vector (e.g. 2D/3D translation for grab brush) */
	float dvec[3];
	
	/* brush geometry (bounding box) */
	rcti brush_rect;
	
	/* Custom data for certain brushes */
	/* - map from bGPDstroke's to structs containing custom data about those strokes */
	GHash *stroke_customdata;
	/* - general customdata */
	void *customdata;
	
	
	/* Timer for in-place accumulation of brush effect */
	wmTimer *timer;
	bool timerTick; /* is this event from a timer */
} tGP_BrushEditData;


/* Callback for performing some brush operation on a single point */
typedef bool (*GP_BrushApplyCb)(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                const int radius, const int co[2]);

/* ************************************************ */
/* Utility Functions */

/* Context ---------------------------------------- */

/* Get the sculpting settings */
static GP_BrushEdit_Settings *gpsculpt_get_settings(Scene *scene)
{
	return &scene->toolsettings->gp_sculpt;
}

/* Get the active brush */
static GP_EditBrush_Data *gpsculpt_get_brush(Scene *scene)
{
	GP_BrushEdit_Settings *gset = &scene->toolsettings->gp_sculpt;
	return &gset->brush[gset->brushtype];
}

/* Brush Operations ------------------------------- */

/* Invert behaviour of brush? */
static bool gp_brush_invert_check(tGP_BrushEditData *gso)
{
	/* The basic setting is the brush's setting (from the panel) */
	bool invert = ((gso->brush->flag & GP_EDITBRUSH_FLAG_INVERT) != 0);
	
	/* During runtime, the user can hold down the Ctrl key to invert the basic behaviour */
	if (gso->flag & GP_EDITBRUSH_FLAG_INVERT) {
		invert ^= true;
	}
		
	return invert;
}

/* Compute strength of effect */
static float gp_brush_influence_calc(tGP_BrushEditData *gso, const int radius, const int co[2])
{
	GP_EditBrush_Data *brush = gso->brush;
	
	/* basic strength factor from brush settings */
	float influence = brush->strength;
	
	/* use pressure? */
	if (brush->flag & GP_EDITBRUSH_FLAG_USE_PRESSURE) {
		influence *= gso->pressure;
	}
	
	/* distance fading */
	if (brush->flag & GP_EDITBRUSH_FLAG_USE_FALLOFF) {
		float distance = (float)len_v2v2_int(gso->mval, co);
		float fac;
		
		CLAMP(distance, 0.0f, (float)radius);
		fac = 1.0f - (distance / (float)radius);
		
		influence *= fac;
	}
	
	/* return influence */
	return influence;
}

/* ************************************************ */
/* Brush Callbacks */
/* This section defines the callbacks used by each brush to perform their magic.
 * These are called on each point within the brush's radius.
 */

/* ----------------------------------------------- */
/* Smooth Brush */

/* A simple (but slower + inaccurate) smooth-brush implementation to test the algorithm for stroke smoothing */
static bool gp_brush_smooth_apply(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                  const int radius, const int co[2])
{
	GP_EditBrush_Data *brush = gso->brush;
	bGPDspoint *pt = &gps->points[i];
	float inf = gp_brush_influence_calc(gso, radius, co);
	float pressure = 0.0f;
	float sco[3] = {0.0f};
	
	/* Do nothing if not enough points to smooth out */
	if (gps->totpoints <= 2) {
		return false;
	}
	
	/* Only affect endpoints by a fraction of the normal strength,
	 * to prevent the stroke from shrinking too much
	 */
	if ((i == 0) || (i == gps->totpoints - 1)) {
		inf *= 0.1f;
	}
	
	/* Compute smoothed coordinate by taking the ones nearby */
	/* XXX: This is potentially slow, and suffers from accumulation error as earlier points are handled before later ones */
	{	
		// XXX: this is hardcoded to look at 2 points on either side of the current one (i.e. 5 items total)
		const int   steps = 2;
		const float average_fac = 1.0f / (float)(steps * 2 + 1);
		int step;
		
		/* add the point itself */
		madd_v3_v3fl(sco, &pt->x, average_fac);
		
		if (brush->flag & GP_EDITBRUSH_FLAG_SMOOTH_PRESSURE) {
			pressure += pt->pressure * average_fac;
		}
		
		/* n-steps before/after current point */
		// XXX: review how the endpoints are treated by this algorithm
		// XXX: falloff measures should also introduce some weighting variations, so that further-out points get less weight
		for (step = 1; step <= steps; step++) {
			bGPDspoint *pt1, *pt2;
			int before = i - step;
			int after  = i + step;
			
			CLAMP_MIN(before, 0);
			CLAMP_MAX(after, gps->totpoints - 1);
			
			pt1 = &gps->points[before];
			pt2 = &gps->points[after];
			
			/* add both these points to the average-sum (s += p[i]/n) */
			madd_v3_v3fl(sco, &pt1->x, average_fac);
			madd_v3_v3fl(sco, &pt2->x, average_fac);
			
			/* do pressure too? */
			if (brush->flag & GP_EDITBRUSH_FLAG_SMOOTH_PRESSURE) {
				pressure += pt1->pressure * average_fac;
				pressure += pt2->pressure * average_fac;
			}
		}
	}
	
	/* Based on influence factor, blend between original and optimal smoothed coordinate */
	interp_v3_v3v3(&pt->x, &pt->x, sco, inf);
	
	if (brush->flag & GP_EDITBRUSH_FLAG_SMOOTH_PRESSURE) {
		pt->pressure = pressure;
	}
	
	return true;
}

/* ----------------------------------------------- */
/* Line Thickness Brush */

/* Make lines thicker or thinner by the specified amounts */
static bool gp_brush_thickness_apply(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                     const int radius, const int co[2])
{
	bGPDspoint *pt = gps->points + i;
	float inf;
	
	/* Compute strength of effect
	 * - We divide the strength by 10, so that users can set "sane" values.
	 *   Otherwise, good default values are in the range of 0.093
	 */
	inf = gp_brush_influence_calc(gso, radius, co) / 10.0f;
	
	/* apply */
	// XXX: this is much too strong, and it should probably do some smoothing with the surrounding stuff
	if (gp_brush_invert_check(gso)) {
		/* make line thinner - reduce stroke pressure */
		pt->pressure -= inf;
	}
	else {
		/* make line thicker - increase stroke pressure */
		pt->pressure += inf;
	}
	
	/* Pressure should stay within [0.0, 1.0]
	 * However, it is nice for volumetric strokes to be able to exceed
	 * the upper end of this range. Therefore, we don't actually clamp
	 * down on the upper end.
	 */
	if (pt->pressure < 0.0f)
		pt->pressure = 0.0f;
	
	return true;
}


/* ----------------------------------------------- */
/* Grab Brush */

/* Custom data per stroke for the Grab Brush
 *
 * This basically defines the strength of the effect for each
 * affected stroke point that was within the initial range of
 * the brush region.
 */
typedef struct tGPSB_Grab_StrokeData {
	/* array of indices to corresponding points in the stroke */
	int   *points;
	/* array of influence weights for each of the included points */
	float *weights;
	
	/* capacity of the arrays */
	int capacity;
	/* actual number of items currently stored */
	int size;
} tGPSB_Grab_StrokeData;

/* initialise custom data for handling this stroke */
static void gp_brush_grab_stroke_init(tGP_BrushEditData *gso, bGPDstroke *gps)
{
	tGPSB_Grab_StrokeData *data = NULL;
	
	BLI_assert(gps->totpoints > 0);
	
	/* Check if there are buffers already (from a prior run) */
	if (BLI_ghash_haskey(gso->stroke_customdata, gps)) {
		/* Ensure that the caches are empty
		 * - Since we reuse these between different strokes, we don't
		 *   want the previous invocation's data polluting the arrays
		 */
		data = BLI_ghash_lookup(gso->stroke_customdata, gps);
		BLI_assert(data != NULL);
		
		data->size = 0; /* minimum requirement - so that we can repopulate again */
		
		memset(data->points, 0, sizeof(int) * data->capacity);
		memset(data->weights, 0, sizeof(float) * data->capacity);
	}
	else {
		/* Create new instance */
		data = MEM_callocN(sizeof(tGPSB_Grab_StrokeData), "GP Stroke Grab Data");
		
		data->capacity = gps->totpoints;
		data->size = 0;
		
		data->points  = MEM_callocN(sizeof(int) * data->capacity, "GP Stroke Grab Indices");
		data->weights = MEM_callocN(sizeof(float) * data->capacity, "GP Stroke Grab Weights");
		
		/* hook up to the cache */
		BLI_ghash_insert(gso->stroke_customdata, gps, data);
	}	
}

/* store references to stroke points in the initial stage */
static bool gp_brush_grab_store_points(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                       const int radius, const int co[2])
{
	tGPSB_Grab_StrokeData *data = BLI_ghash_lookup(gso->stroke_customdata, gps);
	float inf = gp_brush_influence_calc(gso, radius, co);
	
	BLI_assert(data != NULL);
	BLI_assert(data->size < data->capacity);
	
	/* insert this point into the set of affected points */
	data->points[data->size]  = i;
	data->weights[data->size] = inf;
	data->size++;
	
	/* done */
	return true;
}

/* Compute effect vector for grab brush */
static void gp_brush_grab_calc_dvec(tGP_BrushEditData *gso)
{
	/* Convert mouse-movements to movement vector */
	// TODO: incorporate pressure into this?
	// XXX: screen-space strokes in 3D space will suffer!
	if (gso->sa->spacetype == SPACE_VIEW3D) {
		View3D *v3d = gso->sa->spacedata.first;
		RegionView3D *rv3d = gso->ar->regiondata;
		float *rvec = ED_view3d_cursor3d_get(gso->scene, v3d);
		float zfac = ED_view3d_calc_zfac(rv3d, rvec, NULL);
		
		float mval_f[2];
		
		/* convert from 2D screenspace to 3D... */
		mval_f[0] = (float)(gso->mval[0] - gso->mval_prev[0]);
		mval_f[1] = (float)(gso->mval[1] - gso->mval_prev[1]);
		
		ED_view3d_win_to_delta(gso->ar, mval_f, gso->dvec, zfac);
	}
	else {
		/* 2D - just copy */
		// XXX: view2d?
		gso->dvec[0] = (float)(gso->mval[0] - gso->mval_prev[0]);
		gso->dvec[1] = (float)(gso->mval[1] - gso->mval_prev[1]);
		gso->dvec[2] = 0.0f;  /* unused */
	}
}

/* Apply grab transform to all relevant points of the affected strokes */
static void gp_brush_grab_apply_cached(tGP_BrushEditData *gso, bGPDstroke *gps)
{
	tGPSB_Grab_StrokeData *data = BLI_ghash_lookup(gso->stroke_customdata, gps);
	int i;
	
	/* Apply dvec to all of the stored points */
	for (i = 0; i < data->size; i++) {
		bGPDspoint *pt = &gps->points[data->points[i]];
		float delta[3] = {0.0f};
		
		/* adjust the amount of displacement to apply */
		mul_v3_v3fl(delta, gso->dvec, data->weights[i]);
		
		/* apply */
		add_v3_v3(&pt->x, delta);
	}
}

/* free customdata used for handling this stroke */
static void gp_brush_grab_stroke_free(void *ptr)
{
	tGPSB_Grab_StrokeData *data = (tGPSB_Grab_StrokeData *)ptr;
	
	/* free arrays */
	MEM_freeN(data->points);
	MEM_freeN(data->weights);
	
	/* ... and this item itself, since it was also allocated */
	MEM_freeN(data);
}

/* ----------------------------------------------- */
/* Push Brush */
/* NOTE: Depends on gp_brush_grab_calc_dvec() */

static bool gp_brush_push_apply(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                const int radius, const int co[2])
{
	bGPDspoint *pt = gps->points + i;
	float inf = gp_brush_influence_calc(gso, radius, co);
	float delta[3] = {0.0f};
		
	/* adjust the amount of displacement to apply */
	mul_v3_v3fl(delta, gso->dvec, inf);
	
	/* apply */
	add_v3_v3(&pt->x, delta);

	/* done */
	return true;
}

/* ----------------------------------------------- */
/* Pinch Brush */

/* Compute reference midpoint for the brush - this is what we'll be moving towards */
static void gp_brush_calc_midpoint(tGP_BrushEditData *gso)
{
	if (gso->sa->spacetype == SPACE_VIEW3D) {
		/* Convert mouse position to 3D space
		 * See: gpencil_paint.c :: gp_stroke_convertcoords()
		 */
		View3D *v3d = gso->sa->spacedata.first;
		RegionView3D *rv3d = gso->ar->regiondata;
		float *rvec = ED_view3d_cursor3d_get(gso->scene, v3d);
		float zfac = ED_view3d_calc_zfac(rv3d, rvec, NULL);
		
		float mval_f[2] = {UNPACK2(gso->mval)};
		float mval_prj[2];
		float dvec[3];
		
		
		if (ED_view3d_project_float_global(gso->ar, rvec, mval_prj, V3D_PROJ_TEST_NOP) == V3D_PROJ_RET_OK) {
			sub_v2_v2v2(mval_f, mval_prj, mval_f);
			ED_view3d_win_to_delta(gso->ar, mval_f, dvec, zfac);
			sub_v3_v3v3(gso->dvec, rvec, dvec);
		}
		else {
			zero_v3(gso->dvec);
		}
	}
	else {
		/* Just 2D coordinates */
		// XXX: fix View2D offsets later
		gso->dvec[0] = (float)gso->mval[0];
		gso->dvec[1] = (float)gso->mval[1];
		gso->dvec[2] = 0.0f;
	}
}

/* Shrink distance between midpoint and this point... */
static bool gp_brush_pinch_apply(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                 const int radius, const int co[2])
{
	bGPDspoint *pt = gps->points + i;
	float fac, inf;
	float vec[3];
	
	/* Scale down standard influence value to get it more manageable...
	 *  - No damping = Unmanageable at > 0.5 strength
	 *  - Div 10     = Not enough effect
	 *  - Div 5      = Happy medium... (by trial and error)
	 */
	inf = gp_brush_influence_calc(gso, radius, co) / 5.0f;
	
	/* 1) Make this point relative to the cursor/midpoint (dvec) */
	sub_v3_v3v3(vec, &pt->x, gso->dvec);
	
	/* 2) Shrink the distance by pulling the point towards the midpoint
	 *    (0.0 = at midpoint, 1 = at edge of brush region)
	 *                         OR
	 *    Increase the distance (if inverting the brush action!)
	 */
	if (gp_brush_invert_check(gso)) {
		/* Inflate (inverse) */
		fac = 1.0f + (inf * inf); /* squared to temper the effect... */
	}
	else {
		/* Shrink (default) */
		fac = 1.0f - (inf * inf); /* squared to temper the effect... */
	}
	mul_v3_fl(vec, fac);
	
	/* 3) Translate back to original space, with the shrinkage applied */
	add_v3_v3v3(&pt->x, gso->dvec, vec);
	
	/* done */
	return true;
}

/* ----------------------------------------------- */
/* Twist Brush - Rotate Around midpoint */

/* Take the screenspace coordinates of the point, rotate this around the brush midpoint,
 * convert the rotated point and convert it into "data" space
 */

static bool gp_brush_twist_apply(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                 const int radius, const int co[2])
{
	bGPDspoint *pt = gps->points + i;
	float tco[2], rco[2], nco[2];
	float rmat[2][2];
	float angle, inf;
	
	/* Angle to rotate by */
	inf = gp_brush_influence_calc(gso, radius, co);
	angle = DEG2RADF(1.0f) * inf;
	
	if (gp_brush_invert_check(gso)) {
		/* invert angle that we rotate by */
		angle *= -1;
	}
	
	/* Express position of point relative to cursor, ready to rotate */
	tco[0] = (float)(co[0] - gso->mval[0]);
	tco[1] = (float)(co[1] - gso->mval[1]);
	
	/* Rotate point in 2D */
	angle_to_mat2(rmat, angle);
	mul_v2_m2v2(rco, rmat, tco);
	
	/* Convert back to screen-coordinates */
	nco[0] = rco[0] + (float)gso->mval[0];
	nco[1] = rco[1] + (float)gso->mval[1];
	
#if 0
	printf("C: %d %d | P: %d %d -> t: %f %f -> r: %f %f x %f -> %f %f\n",
	       gso->mval[0], gso->mval[1], co[0], co[1],
	       tco[0], tco[1],
	       rco[0], rco[1], angle,
	       nco[0], nco[1]);
#endif
	
	/* convert to dataspace */
	if (gps->flag & GP_STROKE_3DSPACE) {
		/* 3D: Project to 3D space */
		if (gso->sa->spacetype == SPACE_VIEW3D) {
			// XXX: this conversion process sometimes introduces noise to the data -> some parts don't seem to move at all, while others get random offsets
			gp_point_xy_to_3d(&gso->gsc, gso->scene, nco, &pt->x);
		}
		else {
			/* ERROR */
			BLI_assert("3D stroke being sculpted in non-3D view");
		}
	}
	else {
		/* 2D: As-is */
		// XXX: v2d scaling/offset?
		copy_v2_v2(&pt->x, nco);
	}
	
	/* done */
	return true;
}


/* ----------------------------------------------- */
/* Randomize Brush */

/* Apply some random jitter to the point */
static bool gp_brush_randomize_apply(tGP_BrushEditData *gso, bGPDstroke *gps, int i,
                                     const int radius, const int co[2])
{
	bGPDspoint *pt = gps->points + i;
	
	/* Amount of jitter to apply depends on the distance of the point to the cursor,
	 * as well as the strength of the brush
	 */
	const float inf = gp_brush_influence_calc(gso, radius, co) / 2.0f;
	const float fac = BLI_frand() * inf;
	
	/* Jitter is applied perpendicular to the mouse movement vector
	 * - We compute all effects in screenspace (since it's easier)
	 *   and then project these to get the points/distances in
	 *   viewspace as needed
	 */
	float mvec[2], svec[2], nco[2];
	
	/* mouse movement in ints -> floats */
	mvec[0] = (float)(gso->mval[0] - gso->mval_prev[0]);
	mvec[1] = (float)(gso->mval[1] - gso->mval_prev[1]);
	
	/* rotate mvec by 90 degrees... */
	svec[0] = -mvec[1];
	svec[1] =  mvec[0];
	
	//printf("svec = %f %f, ", svec[0], svec[1]);
	
	/* scale the displacement by the random displacement, and apply */
	if (BLI_frand() > 0.5f) {
		mul_v2_fl(svec, -fac);
	}
	else {
		mul_v2_fl(svec, fac);
	}
	
	nco[0] = (float)co[0] + svec[0];
	nco[1] = (float)co[1] + svec[1];
	
	//printf("%f %f (%f), nco = {%f %f}, co = %d %d\n", svec[0], svec[1], fac, nco[0], nco[1], co[0], co[1]);
	
	/* convert to dataspace */
	if (gps->flag & GP_STROKE_3DSPACE) {
		/* 3D: Project to 3D space */
		if (gso->sa->spacetype == SPACE_VIEW3D) {
			gp_point_xy_to_3d(&gso->gsc, gso->scene, nco, &pt->x);
		}
		else {
			/* ERROR */
			BLI_assert("3D stroke being sculpted in non-3D view");
		}
	}
	else {
		/* 2D: As-is */
		// XXX: v2d scaling/offset?
		copy_v2_v2(&pt->x, nco);
	}
	
	/* done */
	return true;
}

/* ************************************************ */
/* Non Callback-Based Brushes */

/* Clone Brush ------------------------------------- */
/* How this brush currently works:
 * - If this is start of the brush stroke, paste immediately under the cursor
 *   by placing the midpoint of the buffer strokes under the cursor now
 *
 * - Otherwise, in:
 *   "Stamp Mode" - Move the newly pasted strokes so that their center follows the cursor
 *   "Continuous" - Repeatedly just paste new copies for where the brush is now
 */

/* Custom state data for clone brush */
typedef struct tGPSB_CloneBrushData {
	/* midpoint of the strokes on the clipboard */
	float buffer_midpoint[3];
	
	/* number of strokes in the paste buffer (and/or to be created each time) */
	size_t totitems;
	
	/* for "stamp" mode, the currently pasted brushes */
	bGPDstroke **new_strokes;
} tGPSB_CloneBrushData;

/* Initialise "clone" brush data */
static void gp_brush_clone_init(bContext *C, tGP_BrushEditData *gso)
{
	tGPSB_CloneBrushData *data;
	bGPDstroke *gps;
	
	/* init custom data */
	gso->customdata = data = MEM_callocN(sizeof(tGPSB_CloneBrushData), "CloneBrushData");
	
	/* compute midpoint of strokes on clipboard */
	for (gps = gp_strokes_copypastebuf.first; gps; gps = gps->next) {
		if (ED_gpencil_stroke_can_use(C, gps)) {
			const float dfac = 1.0f / ((float)gps->totpoints);
			float mid[3] = {0.0f};
			
			bGPDspoint *pt;
			int i;
			
			/* compute midpoint of this stroke */
			for (i = 0, pt = gps->points; i < gps->totpoints; i++, pt++) {
				float co[3];
				
				mul_v3_v3fl(co, &pt->x, dfac);
				add_v3_v3(mid, co);
			}
			
			/* combine this stroke's data with the main data */
			add_v3_v3(data->buffer_midpoint, mid);
			data->totitems++;
		}
	}
	
	/* Divide the midpoint by the number of strokes, to finish averaging it */
	if (data->totitems > 1) {
		mul_v3_fl(data->buffer_midpoint, 1.0f / (float)data->totitems);
	}
	
	/* Create a buffer for storing the current strokes */
	if (1 /*gso->brush->mode == GP_EDITBRUSH_CLONE_MODE_STAMP*/) {
		data->new_strokes = MEM_callocN(sizeof(bGPDstroke *) * data->totitems, "cloned strokes ptr array");
	}
}

/* Free custom data used for "clone" brush */
static void gp_brush_clone_free(tGP_BrushEditData *gso)
{
	tGPSB_CloneBrushData *data = gso->customdata;
	
	/* free strokes array */
	if (data->new_strokes) {
		MEM_freeN(data->new_strokes);
		data->new_strokes = NULL;
	}
	
	/* free the customdata itself */
	MEM_freeN(data);
	gso->customdata = NULL;
}

/* Create new copies of the strokes on the clipboard */
static void gp_brush_clone_add(bContext *C, tGP_BrushEditData *gso)
{
	tGPSB_CloneBrushData *data = gso->customdata;
	
	Scene *scene = gso->scene;
	bGPDlayer *gpl = CTX_data_active_gpencil_layer(C);
	bGPDframe *gpf = gpencil_layer_getframe(gpl, CFRA, true);
	bGPDstroke *gps;
	
	float delta[3];
	size_t strokes_added = 0;
	
	/* Compute amount to offset the points by */
	/* NOTE: This assumes that screenspace strokes are NOT used in the 3D view... */
	
	gp_brush_calc_midpoint(gso); /* this puts the cursor location into gso->dvec */
	sub_v3_v3v3(delta, gso->dvec, data->buffer_midpoint);
	
	/* Copy each stroke into the layer */
	for (gps = gp_strokes_copypastebuf.first; gps; gps = gps->next) {
		if (ED_gpencil_stroke_can_use(C, gps)) {
			bGPDstroke *new_stroke;
			bGPDspoint *pt;
			int i;
			
			/* Make a new stroke */
			new_stroke = MEM_dupallocN(gps);
			
			new_stroke->points = MEM_dupallocN(gps->points);
			new_stroke->next = new_stroke->prev = NULL;
			
			BLI_addtail(&gpf->strokes, new_stroke);
			
			/* Adjust all the stroke's points, so that the strokes
			 * get pasted relative to where the cursor is now
			 */
			for (i = 0, pt = new_stroke->points; i < new_stroke->totpoints; i++, pt++) {
				/* assume that the delta can just be applied, and then everything works */
				add_v3_v3(&pt->x, delta);
			}
			
			/* Store ref for later */
			if ((data->new_strokes) && (strokes_added < data->totitems)) {
				data->new_strokes[strokes_added] = new_stroke;
				strokes_added++;
			}
		}
	}
}

/* Move newly-added strokes around - "Stamp" mode of the Clone brush */
static void gp_brush_clone_adjust(tGP_BrushEditData *gso)
{
	tGPSB_CloneBrushData *data = gso->customdata;
	size_t snum;
	
	/* Compute the amount of movement to apply (overwrites dvec) */
	gp_brush_grab_calc_dvec(gso);
	
	/* For each of the stored strokes, apply the offset to each point */
	/* NOTE: Again this assumes that in the 3D view, we only have 3d space and not screenspace strokes... */
	for (snum = 0; snum < data->totitems; snum++) {
		bGPDstroke *gps = data->new_strokes[snum];
		bGPDspoint *pt;
		int i;
		
		for (i = 0, pt = gps->points; i < gps->totpoints; i++, pt++) {
			if (gso->brush->flag & GP_EDITBRUSH_FLAG_USE_FALLOFF) {
				/* "Smudge" Effect when falloff is enabled */
				float delta[3] = {0.0f};
				int sco[2] = {0};
				float influence;
				
				/* compute influence on point */
				gp_point_to_xy(&gso->gsc, gps, pt, &sco[0], &sco[1]);
				influence = gp_brush_influence_calc(gso, gso->brush->size, sco);
				
				/* adjust the amount of displacement to apply */
				mul_v3_v3fl(delta, gso->dvec, influence);
				
				/* apply */
				add_v3_v3(&pt->x, delta);
			}
			else {
				/* Just apply the offset - All points move perfectly in sync with the cursor */
				add_v3_v3(&pt->x, gso->dvec);
			}
		}
	}
}

/* Entrypoint for applying "clone" brush */
static bool gpsculpt_brush_apply_clone(bContext *C, tGP_BrushEditData *gso)
{
	/* Which "mode" are we operating in? */
	if (gso->first) {
		/* Create initial clones */
		gp_brush_clone_add(C, gso);
	}
	else {
		/* Stamp or Continous Mode */
		if (1 /*gso->brush->mode == GP_EDITBRUSH_CLONE_MODE_STAMP*/) {
			/* Stamp - Proceed to translate the newly added strokes */
			gp_brush_clone_adjust(gso);
		}
		else {
			/* Continuous - Just keep pasting everytime we move */
			/* TODO: The spacing of repeat should be controlled using a "stepsize" or similar property? */
			gp_brush_clone_add(C, gso);
		}
	}
	
	return true;
}

/* ************************************************ */
/* Cursor drawing */

/* Helper callback for drawing the cursor itself */
static void gp_brush_drawcursor(bContext *C, int x, int y, void *UNUSED(customdata))
{
	GP_EditBrush_Data *brush = gpsculpt_get_brush(CTX_data_scene(C));
	
	if (brush) {
		glPushMatrix();
		
		glTranslatef((float)x, (float)y, 0.0f);
		
		/* TODO: toggle between add and remove? */
		glColor4ub(255, 255, 255, 128);
		
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		
		glutil_draw_lined_arc(0.0, M_PI * 2.0, brush->size, 40);
		
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
		
		glPopMatrix();
	}
}

/* Turn brush cursor in on/off */
static void gpencil_toggle_brush_cursor(bContext *C, bool enable)
{
	GP_BrushEdit_Settings *gset = gpsculpt_get_settings(CTX_data_scene(C));
	
	if (gset->paintcursor && !enable) {
		/* clear cursor */
		WM_paint_cursor_end(CTX_wm_manager(C), gset->paintcursor);
		gset->paintcursor = NULL;
	}
	else if (enable) {
		/* enable cursor */
		gset->paintcursor = WM_paint_cursor_activate(CTX_wm_manager(C), 
		                                             NULL, 
		                                             gp_brush_drawcursor, NULL);
	}
}


/* ************************************************ */
/* Header Info for GPencil Sculpt */

static void gpsculpt_brush_header_set(bContext *C, tGP_BrushEditData *gso)
{
	const char *brush_name = NULL;
	char str[256] = "";
	
	RNA_enum_name(rna_enum_gpencil_sculpt_brush_items, gso->brush_type, &brush_name);
	
	BLI_snprintf(str, sizeof(str),
	             IFACE_("GPencil Sculpt: %s Stroke  | LMB to paint | RMB/Escape to Exit"
	                    " | Ctrl to Invert Action | Wheel Up/Down for Size "
	                    " | Shift-Wheel Up/Down for Strength"),
	             (brush_name) ? brush_name : "<?>");
	
	ED_area_headerprint(CTX_wm_area(C), str);
}

/* ************************************************ */
/* Grease Pencil Sculpting Operator */

/* Init/Exit ----------------------------------------------- */

static bool gpsculpt_brush_init(bContext *C, wmOperator *op)
{
	Scene *scene = CTX_data_scene(C);
	tGP_BrushEditData *gso;
	
	/* setup operator data */
	gso = MEM_callocN(sizeof(tGP_BrushEditData), "tGP_BrushEditData");
	op->customdata = gso;
	
	/* store state */
	gso->settings = gpsculpt_get_settings(scene);
	gso->brush = gpsculpt_get_brush(scene);
	
	gso->brush_type = gso->settings->brushtype;
	
	
	gso->is_painting = false;
	gso->first = true;
	
	gso->gpd = ED_gpencil_data_get_active(C);
	gso->cfra = INT_MAX; /* NOTE: So that first stroke will get handled in init_stroke() */
	
	gso->scene = scene;
	
	gso->sa = CTX_wm_area(C);
	gso->ar = CTX_wm_region(C);
	
	/* initialise custom data for brushes */
	switch (gso->brush_type) {
		case GP_EDITBRUSH_TYPE_CLONE:
		{
			bGPDstroke *gps;
			bool found = false;
			
			/* check that there are some usable strokes in the buffer */
			for (gps = gp_strokes_copypastebuf.first; gps; gps = gps->next) {
				if (ED_gpencil_stroke_can_use(C, gps)) {
					found = true;
					break;
				}
			}
			
			if (found == false) {
				/* STOP HERE! Nothing to paste! */
				BKE_report(op->reports, RPT_ERROR, 
					   "Copy some strokes to the clipboard before using the Clone brush to paste copies of them");
					   
				MEM_freeN(gso);
				op->customdata = NULL;
				return false;
			}
			else {
				/* initialise customdata */
				gp_brush_clone_init(C, gso);
			}
			break;
		}
		
		case GP_EDITBRUSH_TYPE_GRAB:
		{
			/* initialise the cache needed for this brush */
			gso->stroke_customdata = BLI_ghash_ptr_new("GP Grab Brush - Strokes Hash");
			break;
		}
			
		/* Others - No customdata needed */
		default:
			break;
	}
	
	
	/* setup space conversions */
	gp_point_conversion_init(C, &gso->gsc);
	
	/* update header */
	gpsculpt_brush_header_set(C, gso);
	
	/* setup cursor drawing */
	WM_cursor_modal_set(CTX_wm_window(C), BC_CROSSCURSOR);
	gpencil_toggle_brush_cursor(C, true);
	
	return true;
}

static void gpsculpt_brush_exit(bContext *C, wmOperator *op)
{
	tGP_BrushEditData *gso = op->customdata;
	wmWindow *win = CTX_wm_window(C);
	
	/* free brush-specific data */
	switch (gso->brush_type) {
		case GP_EDITBRUSH_TYPE_GRAB:
		{
			/* Free per-stroke customdata
			 * - Keys don't need to be freed, as those are the strokes
			 * - Values assigned to those keys do, as they are custom structs
			 */
			BLI_ghash_free(gso->stroke_customdata, NULL, gp_brush_grab_stroke_free);
			break;
		}
		
		case GP_EDITBRUSH_TYPE_CLONE:
		{
			/* Free customdata */
			gp_brush_clone_free(gso);
			break;
		}
		
		default:
			break;
	}
	
	/* unregister timer (only used for realtime) */
	if (gso->timer) {
		WM_event_remove_timer(CTX_wm_manager(C), win, gso->timer);
	}

	/* disable cursor and headerprints */
	ED_area_headerprint(CTX_wm_area(C), NULL);
	WM_cursor_modal_restore(win);
	gpencil_toggle_brush_cursor(C, false);
	
	/* free operator data */
	MEM_freeN(gso);
	op->customdata = NULL;
}

/* poll callback for stroke sculpting operator(s) */
static int gpsculpt_brush_poll(bContext *C)
{
	/* NOTE: this is a bit slower, but is the most accurate... */
	return CTX_DATA_COUNT(C, editable_gpencil_strokes) != 0;
}

/* Init Sculpt Stroke ---------------------------------- */

static void gpsculpt_brush_init_stroke(tGP_BrushEditData *gso)
{
	Scene *scene = gso->scene;
	bGPdata *gpd = gso->gpd;
	bGPDlayer *gpl;
	int cfra = CFRA;
	
	/* only try to add a new frame if this is the first stroke, or the frame has changed */
	if ((gpd == NULL) || (cfra == gso->cfra))
		return;
	
	/* go through each layer, and ensure that we've got a valid frame to use */
	for (gpl = gpd->layers.first; gpl; gpl = gpl->next) {
		/* only editable and visible layers are considered */
		if ((gpl->flag & (GP_LAYER_HIDE | GP_LAYER_LOCKED)) == 0 &&
		    (gpl->actframe != NULL))
		{
			bGPDframe *gpf = gpl->actframe;
			
			/* Make a new frame to work on if the layer's frame and the current scene frame don't match up 
			 * - This is useful when animating as it saves that "uh-oh" moment when you realize you've
			 *   spent too much time editing the wrong frame...
			 */
			// XXX: should this be allowed when framelock is enabled?
			if (gpf->framenum != cfra) {
				gpencil_frame_addcopy(gpl, cfra);
			}
		}
	}
	
	/* save off new current frame, so that next update works fine */
	gso->cfra = cfra;
}

/* Apply ----------------------------------------------- */

/* Apply brush operation to points in this stroke */
static bool gpsculpt_brush_do_stroke(tGP_BrushEditData *gso, bGPDstroke *gps, GP_BrushApplyCb apply)
{
	GP_SpaceConversion *gsc = &gso->gsc;
	rcti *rect = &gso->brush_rect;
	const int radius = gso->brush->size;
	
	bGPDspoint *pt1, *pt2;
	int pc1[2] = {0};
	int pc2[2] = {0};
	int i;
	bool include_last = false;
	bool changed = false;
	
	if (gps->totpoints == 1) {
		gp_point_to_xy(gsc, gps, gps->points, &pc1[0], &pc1[1]);
		
		/* do boundbox check first */
		if ((!ELEM(V2D_IS_CLIPPED, pc1[0], pc1[1])) && BLI_rcti_isect_pt(rect, pc1[0], pc1[1])) {
			/* only check if point is inside */
			if (len_v2v2_int(gso->mval, pc1) <= radius) {
				/* apply operation to this point */
				changed = apply(gso, gps, 0, radius, pc1);
			}
		}
	}
	else {
		/* Loop over the points in the stroke, checking for intersections 
		 *  - an intersection means that we touched the stroke
		 */
		for (i = 0; (i + 1) < gps->totpoints; i++) {
			/* Get points to work with */
			pt1 = gps->points + i;
			pt2 = gps->points + i + 1;
			
			/* Skip if neither one is selected (and we are only allowed to edit/consider selected points) */
			if (gso->settings->flag & GP_BRUSHEDIT_FLAG_SELECT_MASK) {
				if (!(pt1->flag & GP_SPOINT_SELECT) && !(pt2->flag & GP_SPOINT_SELECT)) {
					include_last = false;
					continue;
				}
			}
			
			gp_point_to_xy(gsc, gps, pt1, &pc1[0], &pc1[1]);
			gp_point_to_xy(gsc, gps, pt2, &pc2[0], &pc2[1]);
			
			/* Check that point segment of the boundbox of the selection stroke */
			if (((!ELEM(V2D_IS_CLIPPED, pc1[0], pc1[1])) && BLI_rcti_isect_pt(rect, pc1[0], pc1[1])) ||
			    ((!ELEM(V2D_IS_CLIPPED, pc2[0], pc2[1])) && BLI_rcti_isect_pt(rect, pc2[0], pc2[1])))
			{
				/* Check if point segment of stroke had anything to do with
				 * eraser region  (either within stroke painted, or on its lines)
				 *  - this assumes that linewidth is irrelevant
				 */
				if (gp_stroke_inside_circle(gso->mval, gso->mval_prev, radius, pc1[0], pc1[1], pc2[0], pc2[1])) {
					/* Apply operation to these points */
					bool ok = false;
					
					/* To each point individually... */
					ok = apply(gso, gps, i, radius, pc1);
					
					/* Only do the second point if this is the last segment,
					 * and it is unlikely that the point will get handled
					 * otherwise. 
					 * 
					 * NOTE: There is a small risk here that the second point wasn't really
					 *       actually in-range. In that case, it only got in because
					 *       the line linking the points was!
					 */
					if (i + 1 == gps->totpoints - 1) {
						ok |= apply(gso, gps, i + 1, radius, pc2);
						include_last = false;
					}
					else {
						include_last = true;
					}
					
					changed |= ok;
				}
				else if (include_last) {
					/* This case is for cases where for whatever reason the second vert (1st here) doesn't get included
					 * because the whole edge isn't in bounds, but it would've qualified since it did with the
					 * previous step (but wasn't added then, to avoid double-ups) 
					 */
					changed |= apply(gso, gps, i, radius, pc1);
					include_last = false;
				}
			}
		}
	}
	
	return changed;
}

/* Perform two-pass brushes which modify the existing strokes */
static bool gpsculpt_brush_apply_standard(bContext *C, tGP_BrushEditData *gso)
{
	bool changed = false;
	
	/* Calculate brush-specific data which applies equally to all points */
	switch (gso->brush_type) {
		case GP_EDITBRUSH_TYPE_GRAB: /* Grab points */
		case GP_EDITBRUSH_TYPE_PUSH: /* Push points */
		{
			/* calculate amount of displacement to apply */
			gp_brush_grab_calc_dvec(gso);
			break;
		}
		
		case GP_EDITBRUSH_TYPE_PINCH: /* Pinch points */
		//case GP_EDITBRUSH_TYPE_TWIST: /* Twist points around midpoint */
		{
			/* calculate midpoint of the brush (in data space) */
			gp_brush_calc_midpoint(gso);
			break;
		}
		
		case GP_EDITBRUSH_TYPE_RANDOMIZE: /* Random jitter */
		{
			/* compute the displacement vector for the cursor (in data space) */
			gp_brush_grab_calc_dvec(gso);
			break;
		}
		
		default:
			break;
	}
	
	
	/* Find visible strokes, and perform operations on those if hit */
	CTX_DATA_BEGIN(C, bGPDstroke *, gps, editable_gpencil_strokes)
	{
		switch (gso->brush_type) {
			case GP_EDITBRUSH_TYPE_SMOOTH: /* Smooth strokes */
			{
				changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_smooth_apply);
				break;
			}
			
			case GP_EDITBRUSH_TYPE_THICKNESS: /* Adjust stroke thickness */
			{
				changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_thickness_apply);
				break;
			}
			
			case GP_EDITBRUSH_TYPE_GRAB: /* Grab points */
			{
				if (gso->first) {
					/* First time this brush stroke is being applied:
					 * 1) Prepare data buffers (init/clear) for this stroke
					 * 2) Use the points now under the cursor
					 */
					gp_brush_grab_stroke_init(gso, gps);
					changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_grab_store_points);
				}
				else {
					/* Apply effect to the stored points */
					gp_brush_grab_apply_cached(gso, gps);
					changed |= true;
				}
				break;
			}
			
			case GP_EDITBRUSH_TYPE_PUSH: /* Push points */
			{
				changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_push_apply);
				break;
			}
			
			case GP_EDITBRUSH_TYPE_PINCH: /* Pinch points */
			{
				changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_pinch_apply);
				break;
			}
			
			case GP_EDITBRUSH_TYPE_TWIST: /* Twist points around midpoint */
			{
				changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_twist_apply);
				break;
			}
			
			case GP_EDITBRUSH_TYPE_RANDOMIZE: /* Apply jitter */
			{
				changed |= gpsculpt_brush_do_stroke(gso, gps, gp_brush_randomize_apply);
				break;
			}
			
			default:
				printf("ERROR: Unknown type of GPencil Sculpt brush - %u\n", gso->brush_type);
				break;
		}
	}
	CTX_DATA_END;
	
	return changed;
}

/* Calculate settings for applying brush */
static void gpsculpt_brush_apply(bContext *C, wmOperator *op, PointerRNA *itemptr)
{
	tGP_BrushEditData *gso = op->customdata;
	const int radius = gso->brush->size;
	float mousef[2];
	int mouse[2];
	bool changed = false;
	
	/* Get latest mouse coordinates */
	RNA_float_get_array(itemptr, "mouse", mousef);
	gso->mval[0] = mouse[0] = (int)(mousef[0]);
	gso->mval[1] = mouse[1] = (int)(mousef[1]);
	
	gso->pressure = RNA_float_get(itemptr, "pressure");
	
	if (RNA_boolean_get(itemptr, "pen_flip"))
		gso->flag |= GP_EDITBRUSH_FLAG_INVERT;
	else
		gso->flag &= ~GP_EDITBRUSH_FLAG_INVERT;
	
	
	/* Store coordinates as reference, if operator just started running */
	if (gso->first) {
		gso->mval_prev[0]  = gso->mval[0];
		gso->mval_prev[1]  = gso->mval[1];
		gso->pressure_prev = gso->pressure;
	}
	
	/* Update brush_rect, so that it represents the bounding rectangle of brush */
	gso->brush_rect.xmin = mouse[0] - radius;
	gso->brush_rect.ymin = mouse[1] - radius;
	gso->brush_rect.xmax = mouse[0] + radius;
	gso->brush_rect.ymax = mouse[1] + radius;
	
	
	/* Apply brush */
	if (gso->brush_type == GP_EDITBRUSH_TYPE_CLONE) {
		changed = gpsculpt_brush_apply_clone(C, gso);
	}
	else {
		changed = gpsculpt_brush_apply_standard(C, gso);
	}
	
	
	/* Updates */
	if (changed) {
		WM_event_add_notifier(C, NC_GPENCIL | ND_DATA | NA_EDITED, NULL);
	}
	
	/* Store values for next step */
	gso->mval_prev[0]  = gso->mval[0];
	gso->mval_prev[1]  = gso->mval[1];
	gso->pressure_prev = gso->pressure;
	gso->first = false;
}

/* Running --------------------------------------------- */

/* helper - a record stroke, and apply paint event */
static void gpsculpt_brush_apply_event(bContext *C, wmOperator *op, const wmEvent *event)
{
	tGP_BrushEditData *gso = op->customdata;
	PointerRNA itemptr;
	float mouse[2];
	int tablet = 0;
	
	mouse[0] = event->mval[0] + 1;
	mouse[1] = event->mval[1] + 1;
	
	/* fill in stroke */
	RNA_collection_add(op->ptr, "stroke", &itemptr);
	
	RNA_float_set_array(&itemptr, "mouse", mouse);
	RNA_boolean_set(&itemptr, "pen_flip", event->ctrl != false);
	RNA_boolean_set(&itemptr, "is_start", gso->first);
	
	/* handle pressure sensitivity (which is supplied by tablets) */
	if (event->tablet_data) {
		const wmTabletData *wmtab = event->tablet_data;
		float pressure = wmtab->Pressure;
		
		tablet = (wmtab->Active != EVT_TABLET_NONE);
		
		/* special exception here for too high pressure values on first touch in
		 * windows for some tablets: clamp the values to be sane
		 */
		if (tablet && (pressure >= 0.99f)) {
			pressure = 1.0f;
		}
		RNA_float_set(&itemptr, "pressure", pressure);
	}
	else {
		RNA_float_set(&itemptr, "pressure", 1.0f);
	}
	
	/* apply */
	gpsculpt_brush_apply(C, op, &itemptr);
}

/* reapply */
static int gpsculpt_brush_exec(bContext *C, wmOperator *op)
{
	if (!gpsculpt_brush_init(C, op))
		return OPERATOR_CANCELLED;
	
	RNA_BEGIN(op->ptr, itemptr, "stroke") 
	{
		gpsculpt_brush_apply(C, op, &itemptr);
	}
	RNA_END;
	
	gpsculpt_brush_exit(C, op);
	
	return OPERATOR_FINISHED;
}


/* start modal painting */
static int gpsculpt_brush_invoke(bContext *C, wmOperator *op, const wmEvent *event)
{
	tGP_BrushEditData *gso = NULL;
	const bool is_modal = RNA_boolean_get(op->ptr, "wait_for_input");
	bool needs_timer = false;
	float brush_rate = 0.0f;
	
	/* init painting data */
	if (!gpsculpt_brush_init(C, op))
		return OPERATOR_CANCELLED;
	
	gso = op->customdata;
	
	/* initialise type-specific data (used for the entire session) */
	switch (gso->brush_type) {
		/* Brushes requiring timer... */
		case GP_EDITBRUSH_TYPE_THICKNESS:
			brush_rate = 0.01f; // XXX: hardcoded
			needs_timer = true;
			break;
			
		case GP_EDITBRUSH_TYPE_PINCH:
			brush_rate = 0.001f; // XXX: hardcoded
			needs_timer = true;
			break;
		
		case GP_EDITBRUSH_TYPE_TWIST:
			brush_rate = 0.01f; // XXX: hardcoded
			needs_timer = true;
			break;
			
		default:
			break;
	}
	
	/* register timer for increasing influence by hovering over an area */
	if (needs_timer) {
		gso->timer = WM_event_add_timer(CTX_wm_manager(C), CTX_wm_window(C), TIMER, brush_rate);
	}
	
	/* register modal handler */
	WM_event_add_modal_handler(C, op);
	
	/* start drawing immediately? */
	if (is_modal == false) {
		ARegion *ar = CTX_wm_region(C);
		
		/* ensure that we'll have a new frame to draw on */
		gpsculpt_brush_init_stroke(gso);
		
		/* apply first dab... */
		gso->is_painting = true;
		gpsculpt_brush_apply_event(C, op, event);
		
		/* redraw view with feedback */
		ED_region_tag_redraw(ar);
	}
	
	return OPERATOR_RUNNING_MODAL;
}

/* painting - handle events */
static int gpsculpt_brush_modal(bContext *C, wmOperator *op, const wmEvent *event)
{
	tGP_BrushEditData *gso = op->customdata;
	const bool is_modal = RNA_boolean_get(op->ptr, "wait_for_input");
	bool redraw_region = false;
	
	/* The operator can be in 2 states: Painting and Idling */
	if (gso->is_painting) {
		/* Painting  */
		switch (event->type) {
			/* Mouse Move = Apply somewhere else */
			case MOUSEMOVE:
			case INBETWEEN_MOUSEMOVE:
				/* apply brush effect at new position */
				gpsculpt_brush_apply_event(C, op, event);
				
				/* force redraw, so that the cursor will at least be valid */
				redraw_region = true;
				break;
			
			/* Timer Tick - Only if this was our own timer */
			case TIMER:
				if (event->customdata == gso->timer) {
					gso->timerTick = true;
					gpsculpt_brush_apply_event(C, op, event);
					gso->timerTick = false;
				}
				break;
				
			/* Adjust brush settings */
			/* FIXME: Step increments and modifier keys are hardcoded here! */
			case WHEELUPMOUSE:
			case PADPLUSKEY:
				if (event->shift) {
					/* increase strength */
					gso->brush->strength += 0.05f;
					CLAMP_MAX(gso->brush->strength, 1.0f);
				}
				else {
					/* increase brush size */
					gso->brush->size += 3;
					CLAMP_MAX(gso->brush->size, 300);
				}
					
				redraw_region = true;
				break;
			
			case WHEELDOWNMOUSE: 
			case PADMINUS:
				if (event->shift) {
					/* decrease strength */
					gso->brush->strength -= 0.05f;
					CLAMP_MIN(gso->brush->strength, 0.0f);
				}
				else {
					/* decrease brush size */
					gso->brush->size -= 3;
					CLAMP_MIN(gso->brush->size, 1);
				}
					
				redraw_region = true;
				break;
			
			/* Painting mbut release = Stop painting (back to idle) */
			case LEFTMOUSE:
				//BLI_assert(event->val == KM_RELEASE);
				if (is_modal) {
					/* go back to idling... */
					gso->is_painting = false;
				}
				else {
					/* end sculpt session, since we're not modal */
					gso->is_painting = false;
					
					gpsculpt_brush_exit(C, op);
					return OPERATOR_FINISHED;
				}
				break;
				
			/* Abort painting if any of the usual things are tried */
			case MIDDLEMOUSE:
			case RIGHTMOUSE:
			case ESCKEY:
				gpsculpt_brush_exit(C, op);
				return OPERATOR_FINISHED;
		}
	}
	else {
		/* Idling */
		BLI_assert(is_modal == true);
		
		switch (event->type) {
			/* Painting mbut press = Start painting (switch to painting state) */
			case LEFTMOUSE:
				/* do initial "click" apply */
				gso->is_painting = true;
				gso->first = true;
				
				gpsculpt_brush_init_stroke(gso);
				gpsculpt_brush_apply_event(C, op, event);
				break;
				
			/* Exit modal operator, based on the "standard" ops */
			case RIGHTMOUSE:
			case ESCKEY:
				gpsculpt_brush_exit(C, op);
				return OPERATOR_FINISHED;
				
			/* MMB is often used for view manipulations */
			case MIDDLEMOUSE:
				return OPERATOR_PASS_THROUGH;
			
			/* Mouse movements should update the brush cursor - Just redraw the active region */
			case MOUSEMOVE:
			case INBETWEEN_MOUSEMOVE:
				redraw_region = true;
				break;
			
			/* Adjust brush settings */
			/* FIXME: Step increments and modifier keys are hardcoded here! */
			case WHEELUPMOUSE:
			case PADPLUSKEY:
				if (event->shift) {
					/* increase strength */
					gso->brush->strength += 0.05f;
					CLAMP_MAX(gso->brush->strength, 1.0f);
				}
				else {
					/* increase brush size */
					gso->brush->size += 3;
					CLAMP_MAX(gso->brush->size, 300);
				}
					
				redraw_region = true;
				break;
			
			case WHEELDOWNMOUSE: 
			case PADMINUS:
				if (event->shift) {
					/* decrease strength */
					gso->brush->strength -= 0.05f;
					CLAMP_MIN(gso->brush->strength, 0.0f);
				}
				else {
					/* decrease brush size */
					gso->brush->size -= 3;
					CLAMP_MIN(gso->brush->size, 1);
				}
					
				redraw_region = true;
				break;
			
			/* Change Frame - Allowed */
			case LEFTARROWKEY:
			case RIGHTARROWKEY:
			case UPARROWKEY:
			case DOWNARROWKEY:
				return OPERATOR_PASS_THROUGH;
			
			/* Unhandled event */
			default:
				break;
		}
	}
	
	/* Redraw region? */
	if (redraw_region) {
		ARegion *ar = CTX_wm_region(C);
		ED_region_tag_redraw(ar);
	}
	
	return OPERATOR_RUNNING_MODAL;
}


/* Operator --------------------------------------------- */

void GPENCIL_OT_brush_paint(wmOperatorType *ot)
{
	PropertyRNA *prop;
	
	/* identifiers */
	ot->name = "Stroke Sculpt";
	ot->idname = "GPENCIL_OT_brush_paint";
	ot->description = "Apply tweaks to strokes by painting over the strokes"; // XXX
	
	/* api callbacks */
	ot->exec = gpsculpt_brush_exec;
	ot->invoke = gpsculpt_brush_invoke;
	ot->modal = gpsculpt_brush_modal;
	ot->cancel = gpsculpt_brush_exit;
	ot->poll = gpsculpt_brush_poll;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO | OPTYPE_BLOCKING;

	/* properties */
	RNA_def_collection_runtime(ot->srna, "stroke", &RNA_OperatorStrokeElement, "Stroke", "");
	
	prop = RNA_def_boolean(ot->srna, "wait_for_input", true, "Wait for Input",
	                       "Enter a mini 'sculpt-mode' if enabled, otherwise, exit after drawing a single stroke");
	RNA_def_property_flag(prop, PROP_SKIP_SAVE);
}

/* ************************************************ */

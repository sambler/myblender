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
 * The Original Code is Copyright (C) 2008 Blender Foundation.
 * All rights reserved.
 *
 * 
 * Contributor(s): Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/uvedit/uvedit_intern.h
 *  \ingroup eduv
 */


#ifndef __UVEDIT_INTERN_H__
#define __UVEDIT_INTERN_H__

struct EditFace;
struct EditMesh;
struct Image;
struct MTFace;
struct Object;
struct Scene;
struct SpaceImage;
struct UvElementMap;
struct wmOperatorType;

/* id can be from 0 to 3 */
#define TF_PIN_MASK(id) (TF_PIN1 << id)
#define TF_SEL_MASK(id) (TF_SEL1 << id)

/* geometric utilities */

void uv_center(float uv[][2], float cent[2], int quad);
float uv_area(float uv[][2], int quad);
void uv_copy_aspect(float uv_orig[][2], float uv[][2], float aspx, float aspy);

/* find nearest */

typedef struct NearestHit {
	struct EditFace *efa;
	struct MTFace *tf;

	int vert, uv;
	int edge, vert2;
} NearestHit;

void uv_find_nearest_vert(struct Scene *scene, struct Image *ima, struct EditMesh *em, float co[2], float penalty[2], struct NearestHit *hit);
void uv_find_nearest_edge(struct Scene *scene, struct Image *ima, struct EditMesh *em, float co[2], struct NearestHit *hit);

/* utility tool functions */

struct UvElement *ED_get_uv_element(struct UvElementMap *map, struct EditFace *efa, int index);
void uvedit_live_unwrap_update(struct SpaceImage *sima, struct Scene *scene, struct Object *obedit);

/* smart stitch */

/* object that stores display data for previewing before accepting stitching */
typedef struct StitchPreviewer {
	/* OpenGL requires different calls for Triangles and Quads.
	 * here we'll store the quads of the mesh */
	float *preview_quads;
	/* ...and here we'll store the triangles*/
	float *preview_tris;
	/* preview data. These will be either the previewed vertices or edges depending on stitch mode settings */
	float *preview_stitchable;
	float *preview_unstitchable;
	/* here we'll store the number of triangles and quads to be drawn */
	unsigned int num_tris;
	unsigned int num_quads;
	unsigned int num_stitchable;
	unsigned int num_unstitchable;

	/* store static island Quads */
	float *static_quads;
	/* ...and here we'll store the triangles*/
	float *static_tris;
	unsigned int num_static_tris;
	unsigned int num_static_quads;
} StitchPreviewer;

StitchPreviewer *uv_get_stitch_previewer(void);

/* operators */

void UV_OT_average_islands_scale(struct wmOperatorType *ot);
void UV_OT_cube_project(struct wmOperatorType *ot);
void UV_OT_cylinder_project(struct wmOperatorType *ot);
void UV_OT_from_view(struct wmOperatorType *ot);
void UV_OT_minimize_stretch(struct wmOperatorType *ot);
void UV_OT_pack_islands(struct wmOperatorType *ot);
void UV_OT_reset(struct wmOperatorType *ot);
void UV_OT_sphere_project(struct wmOperatorType *ot);
void UV_OT_unwrap(struct wmOperatorType *ot);
void UV_OT_stitch(struct wmOperatorType *ot);

#endif /* __UVEDIT_INTERN_H__ */


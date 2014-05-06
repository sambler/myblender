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
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */
#ifndef BKE_AMA_H
#define BKE_AMA_H

/** \file BKE_ama.h
 *  \ingroup bke
 */

struct ArrayModifierData;
struct Object;
struct Scene;
struct MVert;
struct MFace;
struct CustomData;

typedef struct IndexMapEntry {
	/* the new vert index that this old vert index maps to */
	int new;
	/* -1 if this vert isn't merged, otherwise the old vert index it
	 * should be replaced with
	 */
	int merge;
	/* 1 if this vert's first copy is merged with the last copy of its
	 * merge target, otherwise 0
	 */
	short merge_final;
} IndexMapEntry;


/* calculations is in local space of deformed object
 * so we store in latmat transform from path coord inside object
 */
typedef struct {
	float dmin[3], dmax[3], dsize, dloc[3];
	float curvespace[4][4], objectspace[4][4], objectspace3[3][3];
	int no_rot_axis;
} CurveDeform;

typedef struct {
	float loc[3];
	float rot[3];
	float scale[3];
	int seed;
} Temp;

int test_index_face_maxvert(struct MFace *mface, struct CustomData *fdata, int mfindex, int nr, int maxvert);
int calc_mapping(IndexMapEntry *indexMap, int oldIndex, int copyNum);

float length_fitcurve(struct ArrayModifierData *amd, struct Scene *scene);
int length_to_count(float length, const float offset[3]);
float count_to_length(int count, const float offset[3]);
float f_rand_max(float max);

void array_offset(const float max_off[3], float rit[3], int prop, int sign, int seed);
void init_mat_oc(const int start, const int end, int *vet_mc);
void init_offset(const int start, const int end, struct ArrayModifierData *ar);
void create_offset(const int n, const int totmat, struct ArrayModifierData *ar, struct Object *ob);
//void array_to_curve(struct Scene *scene, struct Object *cuOb, float (*vertexCos)[3], int numVerts);
void array_to_curve(struct Scene *scene, struct Object *cuOb, struct Object *target, float *vertexCos, float *vec, float *cent);

#endif

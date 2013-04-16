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
 * The Original Code is Copyright (C) 2010 by Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Joseph Eagar
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/editmesh_bvh.c
 *  \ingroup bke
 */

#include "MEM_guardedalloc.h"

#include "DNA_object_types.h"

#include "BLI_math.h"
#include "BLI_bitmap.h"

#include "BKE_DerivedMesh.h"
#include "BKE_editmesh.h"

#include "BKE_editmesh_bvh.h"  /* own include */


struct BMBVHTree {
	BMEditMesh *em;
	BMesh *bm;
	BVHTree *tree;
	float epsilon;
	float maxdist; /* for nearest point search */
	float uv[2];
	
	/* stuff for topological vert search */
	BMVert *v, *curv;
	GHash *gh;
	float curw, curd;
	float co[3], (*cagecos)[3], (*cos)[3];
	int curtag, flag;
};

struct CageUserData {
	int totvert;
	float (*cagecos)[3];
	BLI_bitmap vert_bitmap;
};

static void cage_mapped_verts_callback(void *userData, int index, const float co[3],
                                       const float UNUSED(no_f[3]), const short UNUSED(no_s[3]))
{
	struct CageUserData *data = userData;

	if ((index >= 0 && index < data->totvert) && (!BLI_BITMAP_GET(data->vert_bitmap, index))) {
		BLI_BITMAP_SET(data->vert_bitmap, index);
		copy_v3_v3(data->cagecos[index], co);
	}
}

BMBVHTree *BKE_bmbvh_new(BMEditMesh *em, int flag, struct Scene *scene)
{
	struct BMLoop *(*looptris)[3] = em->looptris;
	BMBVHTree *tree = MEM_callocN(sizeof(*tree), "BMBVHTree");
	DerivedMesh *cage, *final;
	float cos[3][3], (*cagecos)[3] = NULL;
	int i;
	int tottri;

	/* BKE_editmesh_tessface_calc() must be called already */
	BLI_assert(em->tottri != 0 || em->bm->totface == 0);

	/* cage-flag needs scene */
	BLI_assert(scene || !(flag & BMBVH_USE_CAGE));

	tree->em = em;
	tree->bm = em->bm;
	tree->epsilon = FLT_EPSILON * 2.0f;
	tree->flag = flag;

	if (flag & (BMBVH_RESPECT_SELECT)) {
		tottri = 0;
		for (i = 0; i < em->tottri; i++) {
			if (BM_elem_flag_test(looptris[i][0]->f, BM_ELEM_SELECT)) {
				tottri++;
			}
		}
	}
	else if (flag & (BMBVH_RESPECT_HIDDEN)) {
		tottri = 0;
		for (i = 0; i < em->tottri; i++) {
			if (!BM_elem_flag_test(looptris[i][0]->f, BM_ELEM_HIDDEN)) {
				tottri++;
			}
		}
	}
	else {
		tottri = em->tottri;
	}

	tree->tree = BLI_bvhtree_new(tottri, tree->epsilon, 8, 8);
	
	if (flag & BMBVH_USE_CAGE) {
		BLI_bitmap vert_bitmap;
		BMIter iter;
		BMVert *v;
		struct CageUserData data;

		tree->cos = MEM_callocN(sizeof(*tree->cos) * em->bm->totvert, "bmbvh cos");
		BM_ITER_MESH_INDEX (v, &iter, em->bm, BM_VERTS_OF_MESH, i) {
			BM_elem_index_set(v, i); /* set_inline */
			copy_v3_v3(tree->cos[i], v->co);
		}
		em->bm->elem_index_dirty &= ~BM_VERT;


		cage = editbmesh_get_derived_cage_and_final(scene, em->ob, em, &final, CD_MASK_DERIVEDMESH);
		cagecos = MEM_callocN(sizeof(float) * 3 * em->bm->totvert, "bmbvh cagecos");
		
		/* when initializing cage verts, we only want the first cage coordinate for each vertex,
		 * so that e.g. mirror or array use original vertex coordinates and not mirrored or duplicate */
		vert_bitmap = BLI_BITMAP_NEW(em->bm->totvert, __func__);

		data.totvert = em->bm->totvert;
		data.cagecos = cagecos;
		data.vert_bitmap = vert_bitmap;
		
		cage->foreachMappedVert(cage, cage_mapped_verts_callback, &data);

		MEM_freeN(vert_bitmap);
	}
	
	tree->cagecos = cagecos;
	
	for (i = 0; i < em->tottri; i++) {


		if (flag & BMBVH_RESPECT_SELECT) {
			/* note, the arrays wont align now! take care */
			if (!BM_elem_flag_test(em->looptris[i][0]->f, BM_ELEM_SELECT)) {
				continue;
			}
		}
		else if (flag & BMBVH_RESPECT_HIDDEN) {
			/* note, the arrays wont align now! take care */
			if (BM_elem_flag_test(looptris[i][0]->f, BM_ELEM_HIDDEN)) {
				continue;
			}
		}

		if (flag & BMBVH_USE_CAGE) {
			copy_v3_v3(cos[0], cagecos[BM_elem_index_get(looptris[i][0]->v)]);
			copy_v3_v3(cos[1], cagecos[BM_elem_index_get(looptris[i][1]->v)]);
			copy_v3_v3(cos[2], cagecos[BM_elem_index_get(looptris[i][2]->v)]);
		}
		else {
			copy_v3_v3(cos[0], looptris[i][0]->v->co);
			copy_v3_v3(cos[1], looptris[i][1]->v->co);
			copy_v3_v3(cos[2], looptris[i][2]->v->co);
		}

		BLI_bvhtree_insert(tree->tree, i, (float *)cos, 3);
	}
	
	BLI_bvhtree_balance(tree->tree);
	
	return tree;
}

void BKE_bmbvh_free(BMBVHTree *tree)
{
	BLI_bvhtree_free(tree->tree);
	
	if (tree->cagecos)
		MEM_freeN(tree->cagecos);
	if (tree->cos)
		MEM_freeN(tree->cos);
	
	MEM_freeN(tree);
}

/* taken from bvhutils.c */
static void raycallback(void *userdata, int index, const BVHTreeRay *ray, BVHTreeRayHit *hit)
{
	BMBVHTree *tree = userdata;
	BMLoop **ltri = tree->em->looptris[index];
	float dist, uv[2];

	if (isect_ray_tri_v3(ray->origin, ray->direction, ltri[0]->v->co, ltri[1]->v->co, ltri[2]->v->co, &dist, uv) &&
	    (dist < hit->dist))
	{
		hit->dist = dist;
		hit->index = index;

		copy_v3_v3(hit->no, ltri[0]->f->no);

		copy_v3_v3(hit->co, ray->direction);
		normalize_v3(hit->co);
		mul_v3_fl(hit->co, dist);
		add_v3_v3(hit->co, ray->origin);
		
		copy_v2_v2(tree->uv, uv);
	}
}

BMFace *BKE_bmbvh_ray_cast(BMBVHTree *tree, const float co[3], const float dir[3],
                           float *r_dist, float r_hitout[3], float r_cagehit[3])
{
	BVHTreeRayHit hit;
	const float dist = r_dist ? *r_dist : FLT_MAX;

	hit.dist = dist;
	hit.index = -1;

	zero_v2(tree->uv);
	
	BLI_bvhtree_ray_cast(tree->tree, co, dir, 0.0f, &hit, raycallback, tree);
	if (hit.index != -1 && hit.dist != dist) {
		if (r_hitout) {
			if (tree->flag & BMBVH_RETURN_ORIG) {
				BMLoop **ltri = tree->em->looptris[hit.index];
				const float *co1, *co2, *co3;
				int i;

				co1 = ltri[0]->v->co;
				co2 = ltri[1]->v->co;
				co3 = ltri[2]->v->co;
				
				for (i = 0; i < 3; i++) {
					r_hitout[i] = co1[i] + ((co2[i] - co1[i]) * tree->uv[0]) +
					                       ((co3[i] - co1[i]) * tree->uv[1]);
				}
			}
			else {
				copy_v3_v3(r_hitout, hit.co);
			}

			if (r_cagehit) {
				copy_v3_v3(r_cagehit, hit.co);
			}
		}

		if (r_dist) {
			*r_dist = hit.dist;
		}

		return tree->em->looptris[hit.index][0]->f;
	}

	return NULL;
}

BVHTree *BKE_bmbvh_tree_get(BMBVHTree *tree)
{
	return tree->tree;
}

static void vertsearchcallback(void *userdata, int index, const float *UNUSED(co), BVHTreeNearest *hit)
{
	BMBVHTree *tree = userdata;
	BMLoop **ltri = tree->em->looptris[index];
	float dist, maxdist, v[3];
	int i;

	maxdist = tree->maxdist;

	for (i = 0; i < 3; i++) {
		sub_v3_v3v3(v, hit->co, ltri[i]->v->co);

		dist = len_v3(v);
		if (dist < hit->dist && dist < maxdist) {
			copy_v3_v3(hit->co, ltri[i]->v->co);
			copy_v3_v3(hit->no, ltri[i]->v->no);
			hit->dist = dist;
			hit->index = index;
		}
	}
}

BMVert *BKE_bmbvh_find_vert_closest(BMBVHTree *tree, const float co[3], const float maxdist)
{
	BVHTreeNearest hit;

	copy_v3_v3(hit.co, co);
	hit.dist = maxdist * 5;
	hit.index = -1;

	tree->maxdist = maxdist;

	BLI_bvhtree_find_nearest(tree->tree, co, &hit, vertsearchcallback, tree);
	if (hit.dist != FLT_MAX && hit.index != -1) {
		BMLoop **ltri = tree->em->looptris[hit.index];
		float dist, curdist = tree->maxdist;
		int cur = 0, i;

		/* maxdist = tree->maxdist; */  /* UNUSED */

		for (i = 0; i < 3; i++) {
			dist = len_v3v3(hit.co, ltri[i]->v->co);
			if (dist < curdist) {
				cur = i;
				curdist = dist;
			}
		}

		return ltri[cur]->v;
	}

	return NULL;
}

/* UNUSED */
#if 0
static short winding(const float v1[3], const float v2[3], const float v3[3])
/* is v3 to the right of (v1 - v2) ? With exception: v3 == v1 || v3 == v2 */
{
	double inp;

	//inp = (v2[cox] - v1[cox]) * (v1[coy] - v3[coy]) + (v1[coy] - v2[coy]) * (v1[cox] - v3[cox]);
	inp = (v2[0] - v1[0]) * (v1[1] - v3[1]) + (v1[1] - v2[1]) * (v1[0] - v3[0]);

	if (inp < 0.0) {
		return 0;
	}
	else if (inp == 0) {
		if (v1[0] == v3[0] && v1[1] == v3[1]) return 0;
		if (v2[0] == v3[0] && v2[1] == v3[1]) return 0;
	}
	return 1;
}
#endif

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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 by the Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Daniel Dunbar
 *                 Ton Roosendaal,
 *                 Ben Batt,
 *                 Brecht Van Lommel,
 *                 Campbell Barton
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_array.c
 *  \ingroup modifiers
 */


/* Array modifier: duplicates the object multiple times along an axis */

///#include <time.h>
///#include <math.h>
#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_edgehash.h" /// AMA TODO for BLI_edgehash_haskey - still use??
#include "BLI_ghash.h"

#include "DNA_curve_types.h"
#include "DNA_group_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_displist.h"
#include "BKE_curve.h"
#include "BKE_group.h"
#include "BKE_modifier.h"
///#include "BKE_anim.h"

///#include "BLI_rand.h"

#include "BKE_ama.h"
#include "MOD_util.h"

#include "bmesh.h"

#include "depsgraph_private.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/* Due to cyclic dependencies it's possible that curve used for
 * deformation here is not evaluated at the time of evaluating
 * this modifier.
 */
#define CYCLIC_DEPENDENCY_WORKAROUND

static void initData(ModifierData *md)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;

	/* default to 2 duplicates distributed along the x-axis by an
	 * offset of 1 object-width
	 */
	amd->start_cap = amd->mid_cap = amd->end_cap = amd->curve_cap = amd->curve_ob = amd->offset_ob = NULL;
	amd->count = 2;
	zero_v3(amd->offset);
	amd->scale[0] = 1;
	amd->scale[1] = amd->scale[2] = 0;
	amd->length = 0;
	amd->merge_dist = 0.01;
	amd->fit_type = MOD_ARR_FIXEDCOUNT;
	amd->offset_type = MOD_ARR_OFF_RELATIVE;
	amd->flags = 0;
	amd->seed[0] = amd->seed[1] = amd->seed[2] = 0;
	amd->type = MOD_ARR_MOD_NRM;
	amd->mode = !MOD_ARR_MOD_ADV;
	amd->displays = !MOD_ARR_DIS_ADV;
	amd->loc_offset[0] = amd->loc_offset[1] = amd->loc_offset[2] = 0;
	amd->rot_offset[0] = amd->rot_offset[1] = amd->rot_offset[2] = 0;
	amd->scale_offset[0] = amd->scale_offset[1] = amd->scale_offset[2] = 1;
	amd->sign = MOD_ARR_SIGN_P + MOD_ARR_SIGN_L;
	amd->lock = MOD_ARR_LOCK;
	amd->flag_offset = MOD_ARR_PROP + MOD_ARR_LOCAL;
	amd->rays = 1;
	amd->rand_mat = MOD_ARR_MAT;
	amd->mat_ob = MOD_ARR_AR_MAT_RND;
	amd->cont_mat = 1;
	amd->count_mc = 1;
	amd->dist_mc = MOD_ARR_DIST_SEQ;
	amd->rays_dir = MOD_ARR_RAYS_X;
	amd->arr_group = NULL;
	amd->rand_group = !MOD_ARR_RAND_GROUP;
	amd->dist_cu = MOD_ARR_DIST_EVENLY;
	amd->outer_cp = !MOD_ARR_CP_FIRST;
	amd->Mem_Mat_Ob.start_cap = 0;
	amd->Mem_Mat_Ob.end_cap = 0;
}

static void copyData(ModifierData *md, ModifierData *target)
{
#if 0
	ArrayModifierData *amd = (ArrayModifierData *) md;
	ArrayModifierData *tamd = (ArrayModifierData *) target;
#endif
	modifier_copyData_generic(md, target);
}

static void foreachObjectLink(
        ModifierData *md, Object *ob,
        void (*walk)(void *userData, Object *ob, Object **obpoin),
        void *userData)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;

	walk(userData, ob, &amd->start_cap);
	walk(userData, ob, &amd->mid_cap);
	walk(userData, ob, &amd->end_cap);
	walk(userData, ob, &amd->curve_ob);
	walk(userData, ob, &amd->curve_cap);
	walk(userData, ob, &amd->offset_ob);
}

static void updateDepgraph(ModifierData *md, DagForest *forest,
                           struct Scene *UNUSED(scene), Object *UNUSED(ob), DagNode *obNode)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;

	if (amd->start_cap) {
		DagNode *curNode = dag_get_node(forest, amd->start_cap);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->mid_cap) {
		DagNode *curNode = dag_get_node(forest, amd->mid_cap);

		dag_add_relation(forest, curNode, obNode,
				 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->end_cap) {
		DagNode *curNode = dag_get_node(forest, amd->end_cap);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->curve_ob) {
		DagNode *curNode = dag_get_node(forest, amd->curve_ob);
		curNode->eval_flags |= DAG_EVAL_NEED_CURVE_PATH;

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->curve_cap) {
		DagNode *curNode = dag_get_node(forest, amd->curve_cap);

		dag_add_relation(forest, curNode, obNode,
						 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->offset_ob) {
		DagNode *curNode = dag_get_node(forest, amd->offset_ob);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
}

static float vertarray_size(MVert *mvert, int numVerts, int axis)
{
	int i;
	float min_co, max_co;

	/* if there are no vertices, width is 0 */
	if (numVerts == 0) return 0;

	/* find the minimum and maximum coordinates on the desired axis */
	min_co = max_co = mvert->co[axis];
	mvert++;
	for (i = 1; i < numVerts; ++i, ++mvert) {
		if (mvert->co[axis] < min_co) min_co = mvert->co[axis];
		if (mvert->co[axis] > max_co) max_co = mvert->co[axis];
	}

	return max_co - min_co;
}

static int *find_doubles_index_map(BMesh *bm, BMOperator *dupe_op,
                                   const ArrayModifierData *amd,
                                   int *index_map_length)
{
	BMOperator find_op;
	BMOIter oiter;
	BMVert *v, *v2;
	BMElem *ele;
	int *index_map, i;

	BMO_op_initf(bm, &find_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
	             "find_doubles verts=%av dist=%f keep_verts=%s",
	             amd->merge_dist, dupe_op, "geom");

	BMO_op_exec(bm, &find_op);

	i = 0;
	BMO_ITER (ele, &oiter, dupe_op->slots_in, "geom", BM_ALL) {
		BM_elem_index_set(ele, i); /* set_dirty */
		i++;
	}

	BMO_ITER (ele, &oiter, dupe_op->slots_out, "geom.out", BM_ALL) {
		BM_elem_index_set(ele, i); /* set_dirty */
		i++;
	}
	/* above loops over all, so set all to dirty, if this is somehow
	 * setting valid values, this line can be removed - campbell */
	bm->elem_index_dirty |= BM_ALL;

	(*index_map_length) = i;
	index_map = MEM_callocN(sizeof(int) * (*index_map_length), "index_map");

	/*element type argument doesn't do anything here*/
	BMO_ITER (v, &oiter, find_op.slots_out, "targetmap.out", 0) {
		v2 = BMO_iter_map_value_ptr(&oiter);

		index_map[BM_elem_index_get(v)] = BM_elem_index_get(v2) + 1;
	}

	BMO_op_finish(bm, &find_op);

	return index_map;
}

/* Used for start/end cap.
 *
 * this function expects all existing vertices to be tagged,
 * so we can know new verts are not tagged.
 *
 * All verts will be tagged on exit.
 */
static void bm_merge_dm_transform(BMesh *bm, DerivedMesh *dm, float mat[4][4],
                                  const ArrayModifierData *amd,
                                  BMOperator *dupe_op,
                                  BMOpSlot dupe_op_slot_args[BMO_OP_MAX_SLOTS], const char *dupe_slot_name,
                                  BMOperator *weld_op)
{
	const bool is_input = (dupe_op->slots_in == dupe_op_slot_args);
	BMVert *v, *v2, *v3;
	BMIter iter;

	/* Add the DerivedMesh's elements to the BMesh. The pre-existing
	 * elements were already tagged, so the new elements can be
	 * identified by not having the BM_ELEM_TAG flag set. */
	DM_to_bmesh_ex(dm, bm, false);

	if (amd->flags & MOD_ARR_MERGE) {
		/* if merging is enabled, find doubles */
		
		BMOIter oiter;
		BMOperator find_op;
		BMOpSlot *slot_targetmap;

		BMO_op_initf(bm, &find_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
		             is_input ?  /* ugh */
		             "find_doubles verts=%Hv dist=%f keep_verts=%s" :
		             "find_doubles verts=%Hv dist=%f keep_verts=%S",
		             BM_ELEM_TAG, amd->merge_dist,
		             dupe_op, dupe_slot_name);

		/* append the dupe's geom to the findop input verts */
		if (is_input) {
			BMO_slot_buffer_append(&find_op, slots_in, "verts",
			                       dupe_op,  slots_in, dupe_slot_name);
		}
		else if (dupe_op->slots_out == dupe_op_slot_args) {
			BMO_slot_buffer_append(&find_op, slots_in,  "verts",
			                       dupe_op,  slots_out, dupe_slot_name);
		}
		else {
			BLI_assert(0);
		}

		/* transform and tag verts */
		BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
			if (!BM_elem_flag_test(v, BM_ELEM_TAG)) {
				mul_m4_v3(mat, v->co);
				BM_elem_flag_enable(v, BM_ELEM_TAG);
			}
		}

		BMO_op_exec(bm, &find_op);

		slot_targetmap = BMO_slot_get(weld_op->slots_in, "targetmap");

		/* add new merge targets to weld operator */
		BMO_ITER (v, &oiter, find_op.slots_out, "targetmap.out", 0) {
			v2 = BMO_iter_map_value_ptr(&oiter);
			/* check in case the target vertex (v2) is already marked
			 * for merging */
			while ((v3 = BMO_slot_map_elem_get(slot_targetmap, v2))) {
				v2 = v3;
			}
			BMO_slot_map_elem_insert(weld_op, slot_targetmap, v, v2);
		}

		BMO_op_finish(bm, &find_op);
	}
	else {
		/* transform and tag verts */
		BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
			if (!BM_elem_flag_test(v, BM_ELEM_TAG)) {
				mul_m4_v3(mat, v->co);
				BM_elem_flag_enable(v, BM_ELEM_TAG);
			}
		}
	}
}

static void merge_first_last(BMesh *bm,
                             const ArrayModifierData *amd,
                             BMOperator *dupe_first,
                             BMOperator *dupe_last,
                             BMOperator *weld_op)
{
	BMOperator find_op;
	BMOIter oiter;
	BMVert *v, *v2;
	BMOpSlot *slot_targetmap;

	BMO_op_initf(bm, &find_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
	             "find_doubles verts=%s dist=%f keep_verts=%s",
	             dupe_first, "geom", amd->merge_dist,
	             dupe_first, "geom");

	/* append the last dupe's geom to the findop input verts */
	BMO_slot_buffer_append(&find_op,  slots_in,  "verts",
	                       dupe_last, slots_out, "geom.out");

	BMO_op_exec(bm, &find_op);

	/* add new merge targets to weld operator */
	slot_targetmap = BMO_slot_get(weld_op->slots_in, "targetmap");
	BMO_ITER (v, &oiter, find_op.slots_out, "targetmap.out", 0) {
		if (!BMO_slot_map_contains(slot_targetmap, v)) {
			v2 = BMO_iter_map_value_ptr(&oiter);
			BMO_slot_map_elem_insert(weld_op, slot_targetmap, v, v2);
		}
	}

	BMO_op_finish(bm, &find_op);
}

static DerivedMesh *arrayModifier_doArray(ArrayModifierData *amd,
                                          Scene *scene, Object *ob, DerivedMesh *dm,
                                          ModifierApplyFlag flag)
{
	DerivedMesh *result;
	BMesh *bm = DM_to_bmesh(dm, false);
	BMOperator first_dupe_op, dupe_op, old_dupe_op, weld_op;
	BMVert **first_geom = NULL;
	int i, j;
	int index_len = -1;  /* initialize to an invalid value */
	/* offset matrix */
	float offset[4][4];
	float final_offset[4][4];
	float length = amd->length;
	int count = amd->count, maxVerts;
	int *indexMap = NULL;
	DerivedMesh *start_cap = NULL, *end_cap = NULL;
	MVert *src_mvert;
	BMOpSlot *slot_targetmap = NULL;  /* for weld_op */

	/* need to avoid infinite recursion here */
	if (amd->start_cap && amd->start_cap != ob && amd->start_cap->type == OB_MESH)
		start_cap = get_dm_for_modifier(amd->start_cap, flag);
	if (amd->end_cap && amd->end_cap != ob && amd->end_cap->type == OB_MESH)
		end_cap = get_dm_for_modifier(amd->end_cap, flag);

	unit_m4(offset);

	src_mvert = dm->getVertArray(dm);
	maxVerts = dm->getNumVerts(dm);

	if (amd->offset_type & MOD_ARR_OFF_CONST)
		add_v3_v3v3(offset[3], offset[3], amd->offset);
	if (amd->offset_type & MOD_ARR_OFF_RELATIVE) {
		for (j = 0; j < 3; j++)
			offset[3][j] += amd->scale[j] * vertarray_size(src_mvert, maxVerts, j);
	}

	if ((amd->offset_type & MOD_ARR_OFF_OBJ) && (amd->offset_ob)) {
		float obinv[4][4];
		float result_mat[4][4];

		if (ob)
			invert_m4_m4(obinv, ob->obmat);
		else
			unit_m4(obinv);

		mul_serie_m4(result_mat, offset,
					 obinv, amd->offset_ob->obmat,
					 NULL, NULL, NULL, NULL, NULL);
		copy_m4_m4(offset, result_mat);
	}

	if (amd->fit_type == MOD_ARR_FITCURVE && amd->curve_ob) {
		Curve *cu = amd->curve_ob->data;
		if (cu) {
	#ifdef CYCLIC_DEPENDENCY_WORKAROUND
			if (amd->curve_ob->curve_cache == NULL) {
				BKE_displist_make_curveTypes(scene, amd->curve_ob, false);
			}
	#endif

			if (amd->curve_ob->curve_cache && amd->curve_ob->curve_cache->path) {
				float scale = mat4_to_scale(amd->curve_ob->obmat);
				length = scale * amd->curve_ob->curve_cache->path->totdist;
			}
		}
	}

	/* calculate the maximum number of copies which will fit within the
	 * prescribed length */
	if (amd->fit_type == MOD_ARR_FITLENGTH || amd->fit_type == MOD_ARR_FITCURVE) {
		float dist = len_v3(offset[3]);

		if (dist > 1e-6f)
			/* this gives length = first copy start to last copy end
			 * add a tiny offset for floating point rounding errors */
			count = (length + 1e-6f) / dist;
		else
			/* if the offset has no translation, just make one copy */
			count = 1;
	}

	if (count < 1)
		count = 1;

	/* calculate the offset matrix of the final copy (for merging) */
	unit_m4(final_offset);

	for (j = 0; j < count - 1; j++) {
		float tmp_mat[4][4];
		mul_m4_m4m4(tmp_mat, offset, final_offset);
		copy_m4_m4(final_offset, tmp_mat);
	}

	/* BMESH_TODO: bumping up the stack level avoids computing the normals
	 * after every top-level operator execution (and this modifier has the
	 * potential to execute a *lot* of top-level BMOps. There should be a
	 * cleaner way to do this. One possibility: a "mirror" BMOp would
	 * certainly help by compressing it all into one top-level BMOp that
	 * executes a lot of second-level BMOps. */
	BM_mesh_elem_toolflags_ensure(bm);
	BMO_push(bm, NULL);
	bmesh_edit_begin(bm, 0);

	if (amd->flags & MOD_ARR_MERGE) {
		BMO_op_init(bm, &weld_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
					"weld_verts");

		slot_targetmap = BMO_slot_get(weld_op.slots_in, "targetmap");
	}

	BMO_op_initf(bm, &dupe_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
				 "duplicate geom=%avef");
	first_dupe_op = dupe_op;

	for (j = 0; j < count - 1; j++) {
		BMVert *v, *v2, *v3;
		BMOpSlot *geom_slot;
		BMOpSlot *geom_out_slot;
		BMOIter oiter;

		if (j != 0) {
			BMO_op_initf(bm, &dupe_op,
						 (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
						 "duplicate geom=%S", &old_dupe_op, "geom.out");
		}
		BMO_op_exec(bm, &dupe_op);

		geom_slot   = BMO_slot_get(dupe_op.slots_in,  "geom");
		geom_out_slot = BMO_slot_get(dupe_op.slots_out, "geom.out");

		if ((amd->flags & MOD_ARR_MERGEFINAL) && j == 0) {
			int first_geom_bytes = sizeof(BMVert *) * geom_slot->len;

			/* make a copy of the initial geometry ordering so the
			 * last duplicate can be merged into it */
			first_geom = MEM_mallocN(first_geom_bytes, "first_geom");
			memcpy(first_geom, geom_slot->data.buf, first_geom_bytes);
		}

		/* apply transformation matrix */
		BMO_ITER (v, &oiter, dupe_op.slots_out, "geom.out", BM_VERT) {
			mul_m4_v3(offset, v->co);
		}

		if (amd->flags & MOD_ARR_MERGE) {
			/*calculate merge mapping*/
			if (j == 0) {
				indexMap = find_doubles_index_map(bm, &dupe_op,
												  amd, &index_len);
			}

	#define _E(s, i) ((BMVert **)(s)->data.buf)[i]

			/* ensure this is set */
			BLI_assert(index_len != -1);

			for (i = 0; i < index_len; i++) {
				if (!indexMap[i]) continue;

				/* merge v (from 'geom.out') into v2 (from old 'geom') */
				v = _E(geom_out_slot, i - geom_slot->len);
				v2 = _E(geom_slot, indexMap[i] - 1);

				/* check in case the target vertex (v2) is already marked
				 * for merging */
				while ((v3 = BMO_slot_map_elem_get(slot_targetmap, v2))) {
					v2 = v3;
				}

				BMO_slot_map_elem_insert(&weld_op, slot_targetmap, v, v2);
			}

	#undef _E
		}

		/* already copied earlier, but after executation more slot
		 * memory may be allocated */
		if (j == 0)
			first_dupe_op = dupe_op;

		if (j >= 2)
			BMO_op_finish(bm, &old_dupe_op);
		old_dupe_op = dupe_op;
	}

	if ((amd->flags & MOD_ARR_MERGE) &&
		(amd->flags & MOD_ARR_MERGEFINAL) &&
		(count > 1))
	{
		/* Merge first and last copies. Note that we can't use the
		 * indexMap for this because (unless the array is forming a
		 * loop) the offset between first and last is different from
		 * dupe X to dupe X+1. */

		merge_first_last(bm, amd, &first_dupe_op, &dupe_op, &weld_op);
	}

	/* start capping */
	if (start_cap || end_cap) {
		BM_mesh_elem_hflag_enable_all(bm, BM_VERT, BM_ELEM_TAG, false);

		if (start_cap) {
			float startoffset[4][4];
			invert_m4_m4(startoffset, offset);
			bm_merge_dm_transform(bm, start_cap, startoffset, amd,
								  &first_dupe_op, first_dupe_op.slots_in, "geom", &weld_op);
		}

		if (end_cap) {
			float endoffset[4][4];
			mul_m4_m4m4(endoffset, offset, final_offset);
			bm_merge_dm_transform(bm, end_cap, endoffset, amd,
								  &dupe_op, (count == 1) ? dupe_op.slots_in : dupe_op.slots_out,
								  (count == 1) ? "geom" : "geom.out", &weld_op);
		}
	}
	/* done capping */

	/* free remaining dupe operators */
	BMO_op_finish(bm, &first_dupe_op);
	if (count > 2)
		BMO_op_finish(bm, &dupe_op);

	/* run merge operator */
	if (amd->flags & MOD_ARR_MERGE) {
		BMO_op_exec(bm, &weld_op);
		BMO_op_finish(bm, &weld_op);
	}

	/* Bump the stack level back down to match the adjustment up above */
	BMO_pop(bm);

	result = CDDM_from_bmesh(bm, false);

	if ((dm->dirty & DM_DIRTY_NORMALS) ||
		((amd->offset_type & MOD_ARR_OFF_OBJ) && (amd->offset_ob)))
	{
		/* Update normals in case offset object has rotation. */
		result->dirty |= DM_DIRTY_NORMALS;
	}

	BM_mesh_free(bm);

	if (indexMap)
		MEM_freeN(indexMap);
	if (first_geom)
		MEM_freeN(first_geom);

	return result;
	}

static void group_arrayduplilist(ListBase *lb, Scene *scene, Object *ob, int level, int animated)
{
	DupliObject *dob;
	Group *group;
	GroupObject *go;
	float mat[4][4], tmat[4][4], offset[4][4], rot[4][4];
	ModifierData *md;
	int i, cont_rnd;
	float d_alp, alpha=0;

	if(ob->dup_group==NULL) return;
	group= ob->dup_group;

	/* simple preventing of too deep nested groups */
	if(level>MAX_DUPLI_RECUR) return;

	/* handles animated groups, and */
	/* we need to check update for objects that are not in scene... */
	BKE_group_handle_recalc_and_update(NULL, scene, ob, group); /* AMA TODO get ctx instead on NULL */
	animated= animated || BKE_group_is_animated(ob, group);

	for(md=ob->modifiers.first; md; md=md->next) {
		if(md->type == eModifierType_Array) {
			if (md->mode&eModifierMode_Realtime || md->mode&eModifierMode_Render){
				ArrayModifierData *amd = (ArrayModifierData*) md;

				d_alp=0;
				if (amd->rays > 1)
					alpha = (float)6.2831 / amd->rays;
				copy_m4_m4(offset, amd->delta);

				for (i = 0; i < amd->count-1; i++) {
					cont_rnd = 0;
					for(go= group->gobject.first; go; go= go->next) {
						cont_rnd++;
						/* note, if you check on layer here, render goes wrong... it still deforms verts and uses parent imat */
						if(go->ob!=ob) {
							if (amd->rand_group & MOD_ARR_RAND_GROUP) {
								if ((amd->Mem_Ob[i].rand_group_obj != 0) && (amd->Mem_Ob[i].rand_group_obj != cont_rnd))
									continue;
							} else if (amd->rand_group & MOD_ARR_RAND_MAT_GROUP) {
								//assign_material(go->ob, *ob->mat, go->ob->totcol+1);
								//go->ob->actcol = BLI_rand() % ob->totcol;

								//printf("Totcol=%d\n", go->ob->totcol);
								//printf("actcol=%d\n", go->ob->actcol);
							}
							/* Group Dupli Offset, should apply after everything else */
							if (group->dupli_ofs[0] || group->dupli_ofs[1] || group->dupli_ofs[2]) {
								copy_m4_m4(tmat, go->ob->obmat);
								sub_v3_v3v3(tmat[3], tmat[3], group->dupli_ofs);
								mul_m4_m4m4(mat, ob->obmat, tmat);
							} else {
								mul_m4_m4m4(mat, ob->obmat, go->ob->obmat);
							}

							copy_m4_m4(tmat, mat);
							if (amd->rays>1) {
								unit_m4(rot);
								if (amd->rays_dir == MOD_ARR_RAYS_X)
									rotate_m4(rot,'X',d_alp);
								else if (amd->rays_dir == MOD_ARR_RAYS_Y)
									rotate_m4(rot,'Y',d_alp);
								else
									rotate_m4(rot,'Z',d_alp);
								if (d_alp == 0){
									mul_m4_m4m4(mat, tmat, offset);

									copy_m4_m4(tmat, mat);
									mul_m4_m4m4(mat, tmat, rot);
								}
								else{
									mul_m4_m4m4(mat, tmat, offset);

									copy_m4_m4(tmat, mat);
									mul_m4_m4m4(mat, tmat, rot);
								}
							}
							else {
								mul_m4_m4m4(mat, tmat, offset);
							}
							/* Noise */
							if (amd->mode & MOD_ARR_MOD_ADV) {
								if (amd->Mem_Ob[i].transform == 1) {
									float app[4][4];
									unit_m4(app);
									loc_eul_size_to_mat4(app, amd->Mem_Ob[i].loc, amd->Mem_Ob[i].rot, amd->Mem_Ob[i].scale);

									copy_m4_m4(tmat, mat);
									mul_m4_m4m4(mat, tmat, app);
								}
							}

							// dob = new_dupli_object(lb, go->ob, mat, ob->lay, 0, OB_DUPLIARRAY, animated); /// AMA TODO

							if (!(md->mode&eModifierMode_Render))
								dob->no_render = 1;
							else
								dob->no_render = 0;
							/* check the group instance and object layers match, also that the object visible flags are ok. */
							//(G.rendering==0 && dob->ob->restrictflag & OB_RESTRICT_VIEW) || /* AMA TODO from below */
							//(G.rendering && dob->ob->restrictflag & OB_RESTRICT_RENDER)) ||
							if((dob->origlay & group->layer)==0 ||
								!(md->mode&eModifierMode_Realtime) || (!(md->mode&eModifierMode_Editmode) && (ob->mode == OB_MODE_EDIT))) {
								dob->no_draw= 1;
							}
							else {
								dob->no_draw= 0;
							}
							if(go->ob->transflag & OB_DUPLI) {
								copy_m4_m4(dob->ob->obmat, dob->mat);
								/// AMA TODO object_duplilist_recursive(&group->id, scene, go->ob, lb, ob->obmat, level+1, animated);
								///copy_m4_m4(dob->ob->obmat, dob->omat);
							}
						}
					}
					/* Increment for rays */
					if (amd->rays>1) {
						d_alp = d_alp + alpha;
						if (d_alp>6.2831)
							d_alp=0;
					}
					/* Offset for clone group */
					if (d_alp == 0)
						mul_m4_m4m4(offset, amd->delta, offset);
				}
			}
		}
	}
}

static DerivedMesh *applyModifier(ModifierData *md, Object *ob,
                                  DerivedMesh *dm,
                                  ModifierApplyFlag flag)
{
	DerivedMesh *result;
	ArrayModifierData *amd = (ArrayModifierData *) md;

	result = arrayModifier_doArray(amd, md->scene, ob, dm, flag);
	if(result != dm) {
		if(amd->arr_group!=NULL) {
			if (amd->mode & MOD_ARR_MOD_ADV_CLONE) {
				ob->transflag = OB_DUPLIARRAY;
				ob->dup_group = amd->arr_group;
			}
			else {
				ob->transflag = 0;
				ob->dup_group = NULL;
			}
		}
		else {
			ob->transflag = 0;
			ob->dup_group = NULL;
		}
	}
	return result;
}

static void freeData(ModifierData *md)
{
	ArrayModifierData *amd = (ArrayModifierData*) md;

	if (amd)
	{
		if (amd->Mem_Ob)
			MEM_freeN(amd->Mem_Ob);
	}
}

ModifierTypeInfo modifierType_Array = {
	/* name */              "Array",
	/* structName */        "ArrayModifierData",
	/* structSize */        sizeof(ArrayModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsMapping |
	                        eModifierTypeFlag_SupportsEditmode |
	                        eModifierTypeFlag_EnableInEditmode |
	                        eModifierTypeFlag_AcceptsCVs,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          freeData,
	/* isDisabled */        NULL,
	/* updateDepgraph */    updateDepgraph,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};

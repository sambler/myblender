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

/** \file blender/modifiers/intern/MOD_bevel.c
 *  \ingroup modifiers
 */

#include "BLI_utildefines.h"
#include "BLI_math.h"
#include "BLI_string.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_modifier.h"
#include "BKE_mesh.h"
#include "BKE_bmesh.h" /* only for defines */

#include "bmesh.h"

#include "DNA_object_types.h"

#include "MEM_guardedalloc.h"


static void initData(ModifierData *md)
{
	BevelModifierData *bmd = (BevelModifierData *) md;

	bmd->value = 0.1f;
	bmd->res = 1;
	bmd->flags = 0;
	bmd->val_flags = 0;
	bmd->lim_flags = 0;
	bmd->e_flags = 0;
	bmd->bevel_angle = 30;
	bmd->defgrp_name[0] = '\0';
}

static void copyData(ModifierData *md, ModifierData *target)
{
	BevelModifierData *bmd = (BevelModifierData *) md;
	BevelModifierData *tbmd = (BevelModifierData *) target;

	tbmd->value = bmd->value;
	tbmd->res = bmd->res;
	tbmd->flags = bmd->flags;
	tbmd->val_flags = bmd->val_flags;
	tbmd->lim_flags = bmd->lim_flags;
	tbmd->e_flags = bmd->e_flags;
	tbmd->bevel_angle = bmd->bevel_angle;
	BLI_strncpy(tbmd->defgrp_name, bmd->defgrp_name, sizeof(tbmd->defgrp_name));
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *md)
{
	BevelModifierData *bmd = (BevelModifierData *)md;
	CustomDataMask dataMask = 0;

	/* ask for vertexgroups if we need them */
	if (bmd->defgrp_name[0]) dataMask |= CD_MASK_MDEFORMVERT;

	return dataMask;
}

// #define USE_BM_BEVEL_OP_AS_MOD

#ifdef USE_BM_BEVEL_OP_AS_MOD

#define GEOM_MARK   1

/* BMESH_TODO
 *
 * this bevel calls the operator which is missing many of the options
 * which the bevel modifier in trunk has.
 * - no vertex bevel
 * - no weight bevel
 *
 * These will need to be added to the bmesh operator.
 *       - campbell
 *
 * note: this code is very close to MOD_edgesplit.c.
 * note: if 0'd code from trunk included below.
 */
static DerivedMesh *applyModifier(ModifierData *md, struct Object *UNUSED(ob),
                                  DerivedMesh *dm,
                                  ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *result;
	BMesh *bm;
	BMIter iter;
	BMEdge *e;
	BevelModifierData *bmd = (BevelModifierData *) md;
	float threshold = cos((bmd->bevel_angle + 0.00001f) * (float)M_PI / 180.0f);
	const int segments = 16;  /* XXX */

	bm = DM_to_bmesh(dm);

	BM_mesh_normals_update(bm, FALSE);
	BMO_push(bm, NULL);

	if (bmd->lim_flags & BME_BEVEL_ANGLE) {
		BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
			/* check for 1 edge having 2 face users */
			BMLoop *l1, *l2;
			if ((l1 = e->l) &&
			    (l2 = e->l->radial_next) != l1)
			{
				if (dot_v3v3(l1->f->no, l2->f->no) < threshold) {
					BMO_elem_flag_enable(bm, e, GEOM_MARK);
					BMO_elem_flag_enable(bm, e->v1, GEOM_MARK);
					BMO_elem_flag_enable(bm, e->v2, GEOM_MARK);
				}
			}
		}
	}
	else {
		/* crummy, is there a way just to operator on all? - campbell */
		BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
			BMO_elem_flag_enable(bm, e, GEOM_MARK);
			BMO_elem_flag_enable(bm, e->v1, GEOM_MARK);
			BMO_elem_flag_enable(bm, e->v2, GEOM_MARK);
		}
	}

	BMO_op_callf(bm, BMO_FLAG_DEFAULTS,
	             "bevel geom=%fve offset=%f segments=%i",
	             GEOM_MARK, bmd->value, segments);
	BMO_pop(bm);

	result = CDDM_from_bmesh(bm, TRUE);
	BM_mesh_free(bm);

	return result;
}


#else /* from trunk, see note above */

static DerivedMesh *applyModifier(ModifierData *md, Object *UNUSED(ob),
                                  DerivedMesh *derivedData,
                                  ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *result;
	BMesh *bm;

	/*bDeformGroup *def;*/
	int /*i,*/ options, defgrp_index = -1;
	BevelModifierData *bmd = (BevelModifierData *) md;

	options = bmd->flags | bmd->val_flags | bmd->lim_flags | bmd->e_flags;

#if 0
	if ((options & BME_BEVEL_VWEIGHT) && bmd->defgrp_name[0]) {
		defgrp_index = defgroup_name_index(ob, bmd->defgrp_name);
		if (defgrp_index == -1) {
			options &= ~BME_BEVEL_VWEIGHT;
		}
	}
#endif

	bm = DM_to_bmesh(derivedData);
	BME_bevel(bm, bmd->value, bmd->res, options, defgrp_index, DEG2RADF(bmd->bevel_angle), NULL);
	result = CDDM_from_bmesh(bm, TRUE);
	BM_mesh_free(bm);

	/* until we allow for dirty normal flag, always calc,
	 * note: calculating on the CDDM is faster then the BMesh equivalent */
	CDDM_calc_normals(result);

	return result;
}

#endif

static DerivedMesh *applyModifierEM(ModifierData *md, Object *ob,
                                    struct BMEditMesh *UNUSED(editData),
                                    DerivedMesh *derivedData)
{
	return applyModifier(md, ob, derivedData, MOD_APPLY_USECACHE);
}


ModifierTypeInfo modifierType_Bevel = {
	/* name */              "Bevel",
	/* structName */        "BevelModifierData",
	/* structSize */        sizeof(BevelModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsEditmode |
	                        eModifierTypeFlag_EnableInEditmode,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   applyModifierEM,
	/* initData */          initData,
	/* requiredDataMask */  requiredDataMask,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};

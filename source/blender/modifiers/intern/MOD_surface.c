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

/** \file blender/modifiers/intern/MOD_surface.c
 *  \ingroup modifiers
 */


#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_bvhutils.h"
#include "BKE_library.h"
#include "BKE_mesh.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"

#include "MOD_modifiertypes.h"
#include "MOD_util.h"

#include "MEM_guardedalloc.h"


static void initData(ModifierData *md)
{
	SurfaceModifierData *surmd = (SurfaceModifierData *) md;

	surmd->bvhtree = NULL;
}

static void freeData(ModifierData *md)
{
	SurfaceModifierData *surmd = (SurfaceModifierData *) md;

	if (surmd) {
		if (surmd->bvhtree) {
			free_bvhtree_from_mesh(surmd->bvhtree);
			MEM_SAFE_FREE(surmd->bvhtree);
		}

		if (surmd->mesh) {
			BKE_id_free(NULL, surmd->mesh);
			surmd->mesh = NULL;
		}

		MEM_SAFE_FREE(surmd->x);

		MEM_SAFE_FREE(surmd->v);
	}
}

static bool dependsOnTime(ModifierData *UNUSED(md))
{
	return true;
}

static void deformVerts(
        ModifierData *md, const ModifierEvalContext *ctx,
        Mesh *mesh,
        float (*vertexCos)[3],
        int UNUSED(numVerts))
{
	SurfaceModifierData *surmd = (SurfaceModifierData *) md;
	const int cfra = (int)DEG_get_ctime(ctx->depsgraph);

	if (surmd->mesh) {
		BKE_id_free(NULL, surmd->mesh);
	}

	if (mesh) {
		/* Not possible to use get_mesh() in this case as we'll modify its vertices
		 * and get_mesh() would return 'mesh' directly. */
		BKE_id_copy_ex(
		        NULL, (ID *)mesh, (ID **)&surmd->mesh,
		        LIB_ID_CREATE_NO_MAIN |
		        LIB_ID_CREATE_NO_USER_REFCOUNT |
		        LIB_ID_CREATE_NO_DEG_TAG |
		        LIB_ID_COPY_NO_PREVIEW,
		        false);
	}
	else surmd->mesh = MOD_get_mesh_eval(ctx->object, NULL, NULL, NULL, false, false);

	if (!ctx->object->pd) {
		printf("SurfaceModifier deformVerts: Should not happen!\n");
		return;
	}

	if (surmd->mesh) {
		unsigned int numverts = 0, i = 0;
		int init = 0;
		float *vec;
		MVert *x, *v;

		BKE_mesh_apply_vert_coords(surmd->mesh, vertexCos);
		BKE_mesh_calc_normals(surmd->mesh);

		numverts = surmd->mesh->totvert;

		if (numverts != surmd->numverts ||
		    surmd->x == NULL ||
		    surmd->v == NULL ||
		    cfra != surmd->cfra + 1)
		{
			if (surmd->x) {
				MEM_freeN(surmd->x);
				surmd->x = NULL;
			}
			if (surmd->v) {
				MEM_freeN(surmd->v);
				surmd->v = NULL;
			}

			surmd->x = MEM_calloc_arrayN(numverts, sizeof(MVert), "MVert");
			surmd->v = MEM_calloc_arrayN(numverts, sizeof(MVert), "MVert");

			surmd->numverts = numverts;

			init = 1;
		}

		/* convert to global coordinates and calculate velocity */
		for (i = 0, x = surmd->x, v = surmd->v; i < numverts; i++, x++, v++) {
			vec = surmd->mesh->mvert[i].co;
			mul_m4_v3(ctx->object->obmat, vec);

			if (init)
				v->co[0] = v->co[1] = v->co[2] = 0.0f;
			else
				sub_v3_v3v3(v->co, vec, x->co);

			copy_v3_v3(x->co, vec);
		}

		surmd->cfra = cfra;

		if (surmd->bvhtree)
			free_bvhtree_from_mesh(surmd->bvhtree);
		else
			surmd->bvhtree = MEM_callocN(sizeof(BVHTreeFromMesh), "BVHTreeFromMesh");

		if (surmd->mesh->totpoly)
			BKE_bvhtree_from_mesh_get(surmd->bvhtree, surmd->mesh, BVHTREE_FROM_LOOPTRI, 2);
		else
			BKE_bvhtree_from_mesh_get(surmd->bvhtree, surmd->mesh, BVHTREE_FROM_EDGES, 2);
	}
}


ModifierTypeInfo modifierType_Surface = {
	/* name */              "Surface",
	/* structName */        "SurfaceModifierData",
	/* structSize */        sizeof(SurfaceModifierData),
	/* type */              eModifierTypeType_OnlyDeform,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_AcceptsCVs |
	                        eModifierTypeFlag_NoUserAdd,

	/* copyData */          NULL,

	/* deformVerts_DM */    NULL,
	/* deformMatrices_DM */ NULL,
	/* deformVertsEM_DM */  NULL,
	/* deformMatricesEM_DM*/NULL,
	/* applyModifier_DM */  NULL,
	/* applyModifierEM_DM */NULL,

	/* deformVerts */       deformVerts,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     NULL,
	/* applyModifierEM */   NULL,

	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          freeData,
	/* isDisabled */        NULL,
	/* updateDepsgraph */   NULL,
	/* dependsOnTime */     dependsOnTime,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};

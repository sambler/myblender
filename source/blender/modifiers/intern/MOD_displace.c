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

/** \file blender/modifiers/intern/MOD_displace.c
 *  \ingroup modifiers
 */


#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "BLI_utildefines.h"
#include "BLI_ghash.h"
#include "BLI_math.h"
#include "BLI_rand.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_deform.h"
#include "BKE_library.h"
#include "BKE_mesh_mapping.h"
#include "BKE_modifier.h"
#include "BKE_texture.h"

#include "depsgraph_private.h"
#include "MEM_guardedalloc.h"

#include "MOD_util.h"

#include "RE_shader_ext.h"


/* Utils. */

/* Return a "random" float depending on noise_mode. */
static float mode_rng_get_float(DisplaceModifierData *dmd, RNG *rng) {
	const float std_deviation = dmd->noise_std_deviation;

	switch(dmd->noise_mode) {
		/* In constant mode, just return 1, the randomness is provided by "probability" setting. */
		case MOD_DISP_NOISE_MODE_CONSTANT:
			return (BLI_rng_get_float(rng) < 0.5) ? 0.0 : 1.0;
		case MOD_DISP_NOISE_MODE_LINEAR:
			return BLI_rng_get_float(rng);
		case MOD_DISP_NOISE_MODE_GAUSSIAN:
			/* Note: We center our distribution on 0.5. */
			return BLI_rng_get_float_normal(rng, 0.5f, std_deviation);
		default:
			return 0.5f;
	}
}

/* Using two independent rng's allows to tune probability without affecting displace, and vice-versa. */
typedef struct RandData {
	RNG *data;
	RNG *prob;
} RandData;

static void randdata_init(RandData *rng, const int seed_data, const int seed_prob)
{
	rng->data = BLI_rng_new_srandom(seed_data);
	rng->prob = BLI_rng_new_srandom(seed_prob);
}

static void randdata_free(RandData *rng)
{
	BLI_rng_free(rng->data);
	BLI_rng_free(rng->prob);
}

static void reset_looseparts_cache(DisplaceModifierData *dmd, const bool free)
{
	if (free && dmd->loosepart_vertex_map) {
		BLI_assert(dmd->loosepart_vertex_mem != NULL);
		MEM_freeN(dmd->loosepart_vertex_map);
		MEM_freeN(dmd->loosepart_vertex_mem);
	}
	dmd->loosepart_vertex_map = dmd->loosepart_vertex_mem = NULL;
	dmd->cached_lpmap_totedge = dmd->cached_lpmap_totvert = dmd->loosepart_count = -1;
}

/* Common to Texture and Noise Vertex. */

static void displaceModifier_doVertices(DisplaceModifierData *dmd, DerivedMesh *dm,
										float (*vertexCos)[3], const int numVerts, MDeformVert *dvert,
										const int defgrp_index, const float loc_f[3],
										bool (*vget)(DisplaceModifierData *, const int, void *, float[4]), void *data)
{
	int i;
	MVert *mvert;
	const float strength = dmd->strength;
	const float midlevel = dmd->midlevel;
	const int direction = dmd->direction;
	/* In case we do have an offset object, use its distance de mod object as factor. */
	const float nor_fac = ((dmd->flags & MOD_DISP_USE_OBJECT_OFFSET) && dmd->offset_ob) ? len_v3(loc_f) : 1.0f;

	mvert = CDDM_get_verts(dm);

	for (i = 0; i < numVerts; i++) {
		float val[4];
		float weight = 1.0f;
		int j = 3;

		if (dvert) {
			weight = defvert_find_weight(&dvert[i], defgrp_index);
			if (weight == 0.0f) {
				continue;
			}
		}

		/* Get our values (either from texture, or random...). */
		if (!vget(dmd, i, data, val)) {
			continue;
		}

		add_v4_fl(val, -midlevel);

		if (is_zero_v4(val)) {
			continue;
		}

		mul_v4_fl(val, strength * weight);

		if (direction == MOD_DISP_DIR_RGB_XYZ) {
			while (j--) {
				float delta = val[j] * loc_f[j];
				CLAMP(delta, -10000, 10000);
				vertexCos[i][j] += delta;
			}
		}
		else if (direction == MOD_DISP_DIR_VAL_XYZ) {
			while (j--) {
				float delta = val[3] * loc_f[j];
				CLAMP(delta, -10000, 10000);
				vertexCos[i][j] += delta;
			}
		}
		else {  /* if (direction == MOD_DISP_DIR_NOR) */
			while (j--) {
				float delta = val[3] * (mvert[i].no[j] / 32767.0f) * nor_fac;
				CLAMP(delta, -10000, 10000);
				vertexCos[i][j] += delta;
			}
		}
	}
}

/* Texture displace part. */

static bool value_getter_texture(DisplaceModifierData *dmd, const int i, void *data, float r[4])
{
	float (*tex_co)[3] = (float (*)[3])data;
	TexResult texres;

	if (tex_co) {
		texres.nor = NULL;
		BKE_texture_get_value(dmd->modifier.scene, dmd->texture, tex_co[i], &texres, false);
		r[0] = texres.tr;
		r[1] = texres.tg;
		r[2] = texres.tb;
		r[3] = texres.tin;
	}
	else {
		copy_v4_fl(r, 1.0f);  /* Default to white. */
	}

	return true;
}

static void displaceModifier_doTexture(DisplaceModifierData *dmd, Object *ob, DerivedMesh *dm,
									   float (*vertexCos)[3], const int numVerts,
									   MDeformVert *dvert, const int defgrp_index, const float loc_f[3])
{
	float (*tex_co)[3] = NULL;

	if (dmd->texture) {
		tex_co = MEM_callocN(sizeof(*tex_co) * numVerts, __func__);
		get_texture_coords((MappingInfoModifierData *)dmd, ob, dm, vertexCos, tex_co, numVerts);

		modifier_init_texture(dmd->modifier.scene, dmd->texture);
	}

	displaceModifier_doVertices(dmd, dm, vertexCos, numVerts, dvert, defgrp_index, loc_f,
								value_getter_texture, tex_co);

	if (tex_co) {
		MEM_freeN(tex_co);
	}
}

/* Random Vertex part. */

static bool value_getter_random(DisplaceModifierData *dmd, const int UNUSED(i), void *data, float r[4])
{
	RandData *rng = (RandData *)data;
	int i = 4;

	while (i--) {
		r[i] = mode_rng_get_float(dmd, rng->data);
	}

	/* Note: We do it this way, so that we always consume the same amount of randomness, independently from
	 *       noise_probability value.
	 */
	return (BLI_rng_get_float(rng->prob) > dmd->noise_probability) ? false : true;
}

static void displaceModifier_doNoise_Vertices(DisplaceModifierData *dmd, struct Object *UNUSED(ob),
											  struct DerivedMesh *dm,
											  float (*vertexCos)[3], const int numVerts, const int seed,
											  MDeformVert *dvert, const int defgrp_index, const float loc_f[3])
{
	RandData rng;

	randdata_init(&rng, seed, seed + 1);

	displaceModifier_doVertices(dmd, dm, vertexCos, numVerts, dvert, defgrp_index, loc_f,
								value_getter_random, &rng);

	randdata_free(&rng);
}

/* Common to Loose Parts and Whole Mesh. */

void static noise_displace_mesh_part(DisplaceModifierData *dmd, float (*vertexCos)[3], const int numVerts,
									 MDeformVert *dvert, const int defgrp_index, const float loc_f[3],
									 const float rot_f[3], const float scale_f[4], RandData *rng, MeshElemMap *vmap,
									 const bool rotate_or_scale, const bool null_pivot)
{
	const float midlevel = dmd->midlevel;
	const float strength = dmd->strength;

	const int vcount = vmap ? vmap->count : numVerts;
	const float zero_v[3] = {0.0f, 0.0f, 0.0f};

	/* Always get all randomness, to avoid varying results when changing some settings... */
	const float trans[3] = {
		(mode_rng_get_float(dmd, rng->data) - midlevel) * strength * loc_f[0],
		(mode_rng_get_float(dmd, rng->data) - midlevel) * strength * loc_f[1],
		(mode_rng_get_float(dmd, rng->data) - midlevel) * strength * loc_f[2],
	};

	/* We add some tweaking to rotation, so that it would (with midlevel = 0.5 and strength = 1) map to [-pi, pi]. */
	const float rot_corr = 2.0f * (float)M_PI;
	const float rot[3] = {
		(mode_rng_get_float(dmd, rng->data) - midlevel) * strength * rot_f[0] * rot_corr,
		(mode_rng_get_float(dmd, rng->data) - midlevel) * strength * rot_f[1] * rot_corr,
		(mode_rng_get_float(dmd, rng->data) - midlevel) * strength * rot_f[2] * rot_corr,
	};

	/* We add some tweaking to scaling, so that it would never give negative value and, with  midlevel = 0.5 and
	 * strength = 1, it would map to [0, 2] (being centered on 1). Negative scales could be nice, but only as option.
	 * Note we have a mix of per-axis and uniform randomization here.
	 */
	const float scale_ufac = 1.0f + (mode_rng_get_float(dmd, rng->data) - midlevel) * strength * scale_f[3];
	const float scale[3] = {
		max_ff(0.0f, (mode_rng_get_float(dmd, rng->data) - midlevel) * strength * scale_f[0] + scale_ufac),
		max_ff(0.0f, (mode_rng_get_float(dmd, rng->data) - midlevel) * strength * scale_f[1] + scale_ufac),
		max_ff(0.0f, (mode_rng_get_float(dmd, rng->data) - midlevel) * strength * scale_f[2] + scale_ufac),
	};

	float mat[4][4];
	float pivot[3];
	int i;

	/* Check if this loose part should be affected. */
	if (BLI_rng_get_float(rng->prob) > dmd->noise_probability) {
		return;
	}

	if (rotate_or_scale) {
		/* Need a pivot. */
		if (null_pivot && !dvert) {
			/* Used when we get the whole mesh, if no pivot group given, use object's center. */
			zero_v3(pivot);
		}
		else {
			float min[3], max[3];
			INIT_MINMAX(min, max);

			for (i = 0; i < vcount; i++) {
				const int iv = vmap ? vmap->indices[i] : i;

				/* If a vertex group for pivot was specified, only take vertices from this group into account. */
				if (dvert) {
					float weight = defvert_find_weight(&dvert[iv], defgrp_index);
					if (weight <= 0.0f) {
						continue;
					}
				}

				minmax_v3v3_v3(min, max, vertexCos[iv]);
			}

			/* Pivot point is middle of bounding box of loose part. */
			mid_v3_v3v3(pivot, min, max);
		}

		/* Now prepare transform mat4 for this loose part. */
		loc_eul_size_to_mat4(mat, trans, rot, scale);
		transform_pivot_set_m4(mat, pivot);
	}
	else {
		/* Loose part is affected but only with translation, much simpler. */
		unit_m4(mat);
		copy_v3_v3(mat[3], trans);
	}

	for (i = 0; i < vcount; i++) {
		const int iv = vmap ? vmap->indices[i] : i;
		mul_m4_v3(mat, vertexCos[iv]);
	}
}

/* Loose Parts part. */

static void displaceModifier_doNoise_LooseParts(DisplaceModifierData *dmd, struct Object *UNUSED(ob),
												struct DerivedMesh *dm, float (*vertexCos)[3], const int numVerts,
												const int seed, MDeformVert *dvert, const int defgrp_index,
												const float loc_f[3], const float rot_f[3], const float scale_f[4])
{
	RandData rng;

	/* Test if we will need to calculate center of bounding box for pivot */
	const bool rotate_or_scale = !is_zero_v3(rot_f) || !is_zero_v4(scale_f);
	const int numEdges = dm->getNumEdges(dm);

	MeshElemMap *loosepart_vertex_map;
	int i, loosepart_count;

	randdata_init(&rng, seed, seed + 1);

	/* If nothing cached or if totvert or totedge has changed, then update cached loosepart map. */
	if (dmd->loosepart_vertex_map == NULL ||
		(numVerts != dmd->cached_lpmap_totvert) ||
		(numEdges != dmd->cached_lpmap_totedge))
	{
		/* If exists former cached map, then free. */
		reset_looseparts_cache(dmd, true);

		/* Build loosepart map, and store it within modifier data structure */
		BKE_mesh_loosepart_vertex_map_create((MeshElemMap **)(&dmd->loosepart_vertex_map),
											 &dmd->loosepart_vertex_mem, &dmd->loosepart_count,
											 dm->getEdgeArray(dm), numEdges, numVerts);
		dmd->cached_lpmap_totvert = numVerts;
		dmd->cached_lpmap_totedge = numEdges;
	}

	loosepart_vertex_map = (MeshElemMap *)dmd->loosepart_vertex_map;
	loosepart_count = dmd->loosepart_count;

	for(i = 0; i < loosepart_count; i++) {
		noise_displace_mesh_part(dmd, vertexCos, numVerts, dvert, defgrp_index, loc_f,rot_f, scale_f,
								 &rng, &loosepart_vertex_map[i], rotate_or_scale, false);
	}

	randdata_free(&rng);
}

/* Whole Mesh part. */

static void displaceModifier_doNoise_WholeMesh(DisplaceModifierData *dmd, struct Object *UNUSED(ob),
												struct DerivedMesh *UNUSED(dm), float (*vertexCos)[3],
												const int numVerts, const int seed,
												MDeformVert *dvert, const int defgrp_index,
												const float loc_f[3], const float rot_f[3], const float scale_f[4])
{
	RandData rng;

	/* Test if we will need to calculate center of bounding box for pivot */
	const bool rotate_or_scale = !is_zero_v3(rot_f) || !is_zero_v4(scale_f);

	randdata_init(&rng, seed, seed + 1);
	/* For WHOLEMESH mode, BLI_rng does not seem to be random enough...
	 * When objects have very similar names, patterns are visible.
	 * To fix this, we do some extra seeding.  And since it's just once per object it's OK
	 */
	BLI_rng_srandom(rng.data, BLI_rng_get_int(rng.data) );

	if (BLI_rng_get_float(rng.prob) > dmd->noise_probability) {
		return;
	}

	noise_displace_mesh_part(dmd, vertexCos, numVerts, dvert, defgrp_index, loc_f,rot_f, scale_f,
							 &rng, NULL, rotate_or_scale, true);

	randdata_free(&rng);
}

/* Common part. */

static void displaceModifier_do(DisplaceModifierData *dmd, Object *ob, DerivedMesh *dm,
								float (*vertexCos)[3], int numVerts)
{
	float loc_f[3], rot_f[3], scale_f[4];
	MDeformVert *dvert;
	int defgrp_index;

	if (dmd->strength == 0.0f)
		return;

	/* Get final per-axis factors (can also be from other object, in which case loc/rot/scale settings are ignored).
	 * Inspired from MOD_array.c, except here it is not cumulative.
	 */
	if ((dmd->flags & MOD_DISP_USE_OBJECT_OFFSET) && dmd->offset_ob) {
		float obinv[4][4], mat[4][4];

		invert_m4_m4(obinv, ob->obmat);

		/* We want values in modified object's local space, so that there is no changes when both objetcs are
		 * transformed equally.
		 */
		mul_m4_m4m4(mat, obinv, dmd->offset_ob->obmat);

		copy_v3_v3(loc_f, mat[3]);
		mat4_to_eul(rot_f, mat);
		mat4_to_size(scale_f, mat);
		scale_f[3] = 0.0f;  /* No uniform factor in this case. */
	}
	else {
		copy_v3_v3(loc_f, dmd->translation);
		copy_v3_v3(rot_f, dmd->rotation);
		copy_v4_v4(scale_f, dmd->scale);
	}

	modifier_get_vgroup(ob, dm, dmd->defgrp_name, &dvert, &defgrp_index);

	if (dmd->displace_mode == MOD_DISP_MODE_TEXTURE) {
		displaceModifier_doTexture(dmd, ob, dm, vertexCos, numVerts, dvert, defgrp_index, loc_f);
	}
	else {
		int seed = dmd->noise_seed;

		if (dmd->flags & MOD_DISP_NOISE_SEED_OBJID) {
			/* Modify seed with Object ID string, so that copies of same object have different seed. */
			const int hash = BLI_ghashutil_strhash(ob->id.name);
			seed ^= hash;  /* XOR with hash. Seed is then & 255 at init anyway. */
		}

		switch(dmd->displace_mode) {
			case MOD_DISP_MODE_RAND_VERTEX:
				displaceModifier_doNoise_Vertices(dmd, ob, dm, vertexCos, numVerts, seed, dvert, defgrp_index, loc_f);
				break;
			case MOD_DISP_MODE_RAND_LOOSEPARTS:
				displaceModifier_doNoise_LooseParts(dmd, ob, dm, vertexCos, numVerts, seed, dvert, defgrp_index,
													loc_f, rot_f, scale_f);
				break;
			case MOD_DISP_MODE_RAND_WHOLEMESH:
				displaceModifier_doNoise_WholeMesh(dmd, ob, dm, vertexCos, numVerts, seed, dvert, defgrp_index,
												   loc_f, rot_f, scale_f);
				break;
		}
	}
}

/* Common modifiers API. */

static void initData(ModifierData *md)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *) md;

	/* Top level mode choice between TEXTURE and RANDOM */
	dmd->displace_mode = MOD_DISP_MODE_TEXTURE;

	dmd->midlevel = 0.5f;
	dmd->strength = 1.0f;
	dmd->defgrp_name[0] = '\0';
	dmd->direction = MOD_DISP_DIR_VAL_NOR;

	dmd->flags = 0;
	dmd->offset_ob = NULL;
	copy_v3_fl(dmd->translation, 1.0f);
	copy_v3_fl(dmd->rotation, 0.0f);
	copy_v3_fl(dmd->scale, 0.0f);

	/* Fields related to texture displace. */
	dmd->texture = NULL;
	dmd->map_object = NULL;
	dmd->uvlayer_name[0] = '\0';
	dmd->texmapping = MOD_DISP_MAP_LOCAL;

	/* Fields related to random noise displace. */
	dmd->noise_seed = 1;
	dmd->noise_probability = 1.0f;
	dmd->noise_mode = MOD_DISP_NOISE_MODE_LINEAR;
	dmd->noise_std_deviation = 0.5f;

	reset_looseparts_cache(dmd, false);
}

static void copyData(ModifierData *md, ModifierData *target)
{
	DisplaceModifierData *tdmd = (DisplaceModifierData *) target;

	modifier_copyData_generic(md, target);

	if (tdmd->texture) {
		id_us_plus(&tdmd->texture->id);
	}
	reset_looseparts_cache(tdmd, false);
}

static void freeData(ModifierData *md)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *) md;
	if (dmd->texture) {
		id_us_min(&dmd->texture->id);
	}
	reset_looseparts_cache(dmd, true);
}

static CustomDataMask requiredDataMask(Object *UNUSED(ob), ModifierData *md)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *)md;
	CustomDataMask dataMask = 0;

	/* ask for vertexgroups if we need them */
	if (dmd->defgrp_name[0]) {
		dataMask |= CD_MASK_MDEFORMVERT;
	}

	/* ask for UV coordinates if we need them */
	if (dmd->displace_mode == MOD_DISP_MODE_TEXTURE && dmd->texmapping == MOD_DISP_MAP_UV) {
		dataMask |= CD_MASK_MTFACE;
	}

	return dataMask;
}

static bool dependsOnTime(ModifierData *md)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *)md;

	if (dmd->displace_mode == MOD_DISP_MODE_TEXTURE && dmd->texture) {
		return BKE_texture_dependsOnTime(dmd->texture);
	}

	return false;
}

static bool dependsOnNormals(ModifierData *md)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *)md;
	return (ELEM(dmd->displace_mode, MOD_DISP_MODE_TEXTURE, MOD_DISP_MODE_RAND_VERTEX) &&
			dmd->direction == MOD_DISP_DIR_VAL_NOR);
}

static void foreachObjectLink(ModifierData *md, Object *ob,
                              ObjectWalkFunc walk, void *userData)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *) md;

	walk(userData, ob, &dmd->offset_ob);
}

static void foreachIDLink(ModifierData *md, Object *ob,
                          IDWalkFunc walk, void *userData)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *) md;

	walk(userData, ob, (ID **)&dmd->texture);

	foreachObjectLink(md, ob, (ObjectWalkFunc)walk, userData);
}

static void foreachTexLink(ModifierData *md, Object *ob,
                           TexWalkFunc walk, void *userData)
{
	walk(userData, ob, md, "texture");
}

static bool isDisabled(ModifierData *md, int UNUSED(useRenderParams))
{
	DisplaceModifierData *dmd = (DisplaceModifierData *) md;
	if (dmd->strength == 0.0f) {
		return true;
	}
	else if (dmd->displace_mode != MOD_DISP_MODE_TEXTURE && dmd->noise_probability == 0.0f) {
		return true;
	}
	return false;
}

static void updateDepgraph(ModifierData *md, DagForest *forest,
                           struct Main *UNUSED(bmain),
                           struct Scene *UNUSED(scene),
                           Object *UNUSED(ob),
                           DagNode *obNode)
{
	DisplaceModifierData *dmd = (DisplaceModifierData *) md;

	if (dmd->displace_mode == MOD_DISP_MODE_TEXTURE) {
		if (dmd->map_object && dmd->texmapping == MOD_DISP_MAP_OBJECT) {
			DagNode *curNode = dag_get_node(forest, dmd->map_object);
			dag_add_relation(forest, curNode, obNode, DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Displace Modifier");
		}
	}
	else if (dmd->texmapping == MOD_DISP_MAP_GLOBAL) {
		dag_add_relation(forest, obNode, obNode, DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Displace Modifier");
	}

	if ((dmd->flags & MOD_DISP_USE_OBJECT_OFFSET) && dmd->offset_ob) {
		DagNode *curNode = dag_get_node(forest, dmd->offset_ob);
		dag_add_relation(forest, curNode, obNode, DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Displace Modifier");
	}
}

static void deformVerts(ModifierData *md, Object *ob,
                        DerivedMesh *derivedData,
                        float (*vertexCos)[3],
                        int numVerts,
                        ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *dm = get_cddm(ob, NULL, derivedData, vertexCos, dependsOnNormals(md));

	displaceModifier_do((DisplaceModifierData *)md, ob, dm,
	                    vertexCos, numVerts);

	if (dm != derivedData)
		dm->release(dm);
}

static void deformVertsEM(
        ModifierData *md, Object *ob, struct BMEditMesh *editData,
        DerivedMesh *derivedData, float (*vertexCos)[3], int numVerts)
{
	DerivedMesh *dm = get_cddm(ob, editData, derivedData, vertexCos, dependsOnNormals(md));

	displaceModifier_do((DisplaceModifierData *)md, ob, dm,
	                    vertexCos, numVerts);

	if (dm != derivedData)
		dm->release(dm);
}


ModifierTypeInfo modifierType_Displace = {
	/* name */              "Displace",
	/* structName */        "DisplaceModifierData",
	/* structSize */        sizeof(DisplaceModifierData),
	/* type */              eModifierTypeType_OnlyDeform,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsEditmode,

	/* copyData */          copyData,
	/* deformVerts */       deformVerts,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     deformVertsEM,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     NULL,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  requiredDataMask,
	/* freeData */          freeData,
	/* isDisabled */        isDisabled,
	/* updateDepgraph */    updateDepgraph,
	/* dependsOnTime */     dependsOnTime,
	/* dependsOnNormals */	dependsOnNormals,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     foreachIDLink,
	/* foreachTexLink */    foreachTexLink,
};

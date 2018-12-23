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
 * The Original Code is Copyright (C) 2018, Blender Foundation
 * This is a new part of Blender
 *
 * Contributor(s): Antonio Vazquez
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/gpencil_modifiers/intern/MOD_gpenciltime.c
 *  \ingroup modifiers
 */

#include <stdio.h>
#include <string.h>

#include "BLI_utildefines.h"

#include "DNA_meshdata_types.h"
#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_gpencil_modifier_types.h"

#include "BKE_colortools.h"
#include "BKE_context.h"
#include "BKE_deform.h"
#include "BKE_gpencil.h"
#include "BKE_gpencil_modifier.h"

#include "DEG_depsgraph.h"

#include "MOD_gpencil_util.h"
#include "MOD_gpencil_modifiertypes.h"

static void initData(GpencilModifierData *md)
{
	TimeGpencilModifierData *gpmd = (TimeGpencilModifierData *)md;
	gpmd->layername[0] = '\0';
	gpmd->offset = 1;
	gpmd->frame_scale = 1.0f;
	gpmd->flag |= GP_TIME_KEEP_LOOP;
	gpmd->sfra = 1;
	gpmd->efra = 250;
}

static void copyData(const GpencilModifierData *md, GpencilModifierData *target)
{
	BKE_gpencil_modifier_copyData_generic(md, target);
}

static int remapTime(
        struct GpencilModifierData *md, struct Depsgraph *UNUSED(depsgraph),
        struct Scene *scene, struct Object *UNUSED(ob), struct bGPDlayer *gpl, int cfra)
{
	TimeGpencilModifierData *mmd = (TimeGpencilModifierData *)md;
	const bool custom = mmd->flag &	GP_TIME_CUSTOM_RANGE;
	const bool invgpl = mmd->flag & GP_TIME_INVERT_LAYER;
	const bool invpass = mmd->flag & GP_TIME_INVERT_LAYERPASS;
	int sfra = custom ? mmd->sfra : scene->r.sfra;
	int efra = custom ? mmd->efra : scene->r.efra;
	CLAMP_MIN(sfra, 1);
	CLAMP_MIN(efra, 1);
	const int time_range = efra - sfra + 1;
	int offset = mmd->offset;
	int segments = 0;

	/* omit if filter by layer */
	if (mmd->layername[0] != '\0') {
		if (invgpl == false) {
			if (!STREQ(mmd->layername, gpl->info)) {
				return cfra;
			}
		}
		else {
			if (STREQ(mmd->layername, gpl->info)) {
				return cfra;
			}
		}
	}
	/* verify pass */
	if (mmd->layer_pass > 0) {
		if (invpass == false) {
			if (gpl->pass_index != mmd->layer_pass) {
				return cfra;
			}
		}
		else {
			if (gpl->pass_index == mmd->layer_pass) {
				return cfra;
			}
		}
	}

	/* if fix mode, return predefined frame number */
	if (mmd->mode == GP_TIME_MODE_FIX) {
		return offset;
	}

	/* invert current frame number */
	if (mmd->mode == GP_TIME_MODE_REVERSE) {
		cfra = efra - cfra + sfra;
	}

	/* apply frame scale */
	cfra *= mmd->frame_scale;

	/* verify offset never is greater than frame range */
	if (abs(offset) > time_range) {
		offset = offset - ((offset / time_range) * time_range);
	}

	/* verify not outside range if loop is disabled */
	if ((mmd->flag & GP_TIME_KEEP_LOOP) == 0) {
		if (cfra + offset < sfra) {
			return sfra;
		}
		if (cfra + offset > efra) {
			return efra;
		}
	}

	/* check frames before start */
	if (cfra < sfra) {
		segments = ((cfra + sfra) / time_range);
		cfra = cfra + (segments * time_range);
	}

	/* check frames after end */
	if (cfra > efra) {
		segments = ((cfra - sfra) / time_range);
		cfra = cfra - (segments * time_range);
	}

	if (mmd->flag & GP_TIME_KEEP_LOOP) {
		const int nfra = cfra + offset;

		/* if the sum of the cfra is out scene frame range, recalc */
		if (cfra + offset < sfra) {
			const int delta = abs(sfra - nfra);
			return efra - delta + 1;
		}
		else if (cfra + offset > efra) {
			return nfra - efra + sfra - 1;
		}
	}

	return cfra + offset;
}

GpencilModifierTypeInfo modifierType_Gpencil_Time = {
	/* name */              "Time Offset",
	/* structName */        "TimeGpencilModifierData",
	/* structSize */        sizeof(TimeGpencilModifierData),
	/* type */              eGpencilModifierTypeType_Gpencil,
	/* flags */             eGpencilModifierTypeFlag_NoApply,

	/* copyData */          copyData,

	/* deformStroke */      NULL,
	/* generateStrokes */   NULL,
	/* bakeModifier */      NULL,
	/* remapTime */         remapTime,

	/* initData */          initData,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepsgraph */   NULL,
	/* dependsOnTime */     NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
	/* getDuplicationFactor */ NULL,
};

/*
 * $Id$
 *
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
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * Contributor(s): 
 * - Blender Foundation, 2003-2009
 * - Peter Schlaile <peter [at] schlaile [dot] de> 2005/2006
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/seqeffects/seq_transform.c
 *  \ingroup bke
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_math.h" /* for M_PI */
#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"

#include "BKE_sequencer.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

#include "seq_intern.h"

/* **********************************************************************
   TRANSFORM
   ********************************************************************** */
static void init_transform_effect(Sequence *seq)
{
	TransformVars *transform;

	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = MEM_callocN(sizeof(struct TransformVars), "transformvars");

	transform = (TransformVars *)seq->effectdata;

	transform->ScalexIni = 1.0f;
	transform->ScaleyIni = 1.0f;
	transform->ScalexFin = 1.0f;
	transform->ScalexFin = 1.0f;

	transform->xIni=0.0f;
	transform->xFin=0.0f;
	transform->yIni=0.0f;
	transform->yFin=0.0f;

	transform->rotIni=0.0f;
	transform->rotFin=0.0f;

	transform->interpolation=1;
	transform->percent=1;
	transform->uniform_scale=0;
}

static int num_inputs_transform(void)
{
	return 1;
}

static void free_transform_effect(Sequence *seq)
{
	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = NULL;
}

static void copy_transform_effect(Sequence *dst, Sequence *src)
{
	dst->effectdata = MEM_dupallocN(src->effectdata);
}

static void transform_image(int x, int y, struct ImBuf *ibuf1, struct ImBuf *out,
			float scale_x, float scale_y, float translate_x, float translate_y,
			float rotate, int interpolation)
{
	int xo, yo, xi, yi;
	float xt, yt, xr, yr;
	float s,c;

	xo = x;
	yo = y;

	// Rotate
	s= sin(rotate);
	c= cos(rotate);

	for (yi = 0; yi < yo; yi++) {
		for (xi = 0; xi < xo; xi++) {

			//translate point
			xt = xi-translate_x;
			yt = yi-translate_y;

			//rotate point with center ref
			xr =  c*xt + s*yt;
			yr = -s*xt + c*yt;

			//scale point with center ref
			xt = xr / scale_x;
			yt = yr / scale_y;

			//undo reference center point
			xt += (xo / 2.0f);
			yt += (yo / 2.0f);

			//interpolate
			switch(interpolation) {
			case 0:
				neareast_interpolation(ibuf1,out, xt,yt,xi,yi);
				break;
			case 1:
				bilinear_interpolation(ibuf1,out, xt,yt,xi,yi);
				break;
			case 2:
				bicubic_interpolation(ibuf1,out, xt,yt,xi,yi);
				break;
			}
		}
	}
}

static void do_transform(Scene *scene, Sequence *seq, float UNUSED(facf0), int x, int y,
			struct ImBuf *ibuf1,struct ImBuf *out)
{
	TransformVars *transform = (TransformVars *)seq->effectdata;
	float scale_x, scale_y, translate_x, translate_y, rotate_radians;

	// Scale
	if (transform->uniform_scale) {
		scale_x = scale_y = transform->ScalexIni;
	} else {
		scale_x = transform->ScalexIni;
		scale_y = transform->ScaleyIni;
	}

	// Translate
	if(!transform->percent){
		float rd_s = (scene->r.size/100.0f);

		translate_x = transform->xIni*rd_s+(x/2.0f);
		translate_y = transform->yIni*rd_s+(y/2.0f);
	}else{
		translate_x = x*(transform->xIni/100.0f)+(x/2.0f);
		translate_y = y*(transform->yIni/100.0f)+(y/2.0f);
	}

	// Rotate
	rotate_radians = ((float)M_PI*transform->rotIni)/180.0f;

	transform_image(x,y, ibuf1, out, scale_x, scale_y, translate_x, translate_y, rotate_radians, transform->interpolation);
}


static struct ImBuf * do_transform_effect(
	SeqRenderData context, Sequence *seq,float UNUSED(cfra),
	float facf0, float UNUSED(facf1),
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	do_transform(context.scene, seq, facf0,
		context.rectx, context.recty, ibuf1, out);

	return out;
}

/* setup */
void SeqConfigHandle_transform(struct SeqEffectHandle *hndl)
{
	hndl->init = init_transform_effect;
	hndl->num_inputs = num_inputs_transform;
	hndl->free = free_transform_effect;
	hndl->copy = copy_transform_effect;
	hndl->execute = do_transform_effect;
}


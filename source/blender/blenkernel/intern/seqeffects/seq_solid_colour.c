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

/** \file blender/blenkernel/intern/seqeffects/seq_solid_colour.c
 *  \ingroup bke
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"

#include "BKE_sequencer.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

#include "seq_intern.h"

/* **********************************************************************
   SOLID COLOR
   ********************************************************************** */

static void init_solid_color(Sequence *seq)
{
	SolidColorVars *cv;

	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = MEM_callocN(sizeof(struct SolidColorVars), "solidcolor");

	cv = (SolidColorVars *)seq->effectdata;
	cv->col[0] = cv->col[1] = cv->col[2] = 0.5;
}

static int num_inputs_color(void)
{
	return 0;
}

static void free_solid_color(Sequence *seq)
{
	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = NULL;
}

static void copy_solid_color(Sequence *dst, Sequence *src)
{
	dst->effectdata = MEM_dupallocN(src->effectdata);
}

static int early_out_color(struct Sequence *UNUSED(seq), float UNUSED(facf0), float UNUSED(facf1))
{
	return -1;
}

static struct ImBuf * do_solid_color(
	SeqRenderData context, Sequence *seq, float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	SolidColorVars *cv = (SolidColorVars *)seq->effectdata;

	unsigned char *rect;
	float *rect_float;
	int x; /*= context.rectx;*/ /*UNUSED*/
	int y; /*= context.recty;*/ /*UNUSED*/

	if (out->rect) {
		unsigned char col0[3];
		unsigned char col1[3];

		col0[0] = facf0 * cv->col[0] * 255;
		col0[1] = facf0 * cv->col[1] * 255;
		col0[2] = facf0 * cv->col[2] * 255;

		col1[0] = facf1 * cv->col[0] * 255;
		col1[1] = facf1 * cv->col[1] * 255;
		col1[2] = facf1 * cv->col[2] * 255;

		rect = (unsigned char *)out->rect;

		for(y=0; y<out->y; y++) {
			for(x=0; x<out->x; x++, rect+=4) {
				rect[0]= col0[0];
				rect[1]= col0[1];
				rect[2]= col0[2];
				rect[3]= 255;
			}
			y++;
			if (y<out->y) {
				for(x=0; x<out->x; x++, rect+=4) {
					rect[0]= col1[0];
					rect[1]= col1[1];
					rect[2]= col1[2];
					rect[3]= 255;
				}
			}
		}

	} else if (out->rect_float) {
		float col0[3];
		float col1[3];

		col0[0] = facf0 * cv->col[0];
		col0[1] = facf0 * cv->col[1];
		col0[2] = facf0 * cv->col[2];

		col1[0] = facf1 * cv->col[0];
		col1[1] = facf1 * cv->col[1];
		col1[2] = facf1 * cv->col[2];

		rect_float = out->rect_float;

		for(y=0; y<out->y; y++) {
			for(x=0; x<out->x; x++, rect_float+=4) {
				rect_float[0]= col0[0];
				rect_float[1]= col0[1];
				rect_float[2]= col0[2];
				rect_float[3]= 1.0;
			}
			y++;
			if (y<out->y) {
				for(x=0; x<out->x; x++, rect_float+=4) {
					rect_float[0]= col1[0];
					rect_float[1]= col1[1];
					rect_float[2]= col1[2];
					rect_float[3]= 1.0;
				}
			}
		}
	}
	return out;
}

/* setup */
void SeqConfigHandle_solid_colour(struct SeqEffectHandle *hndl)
{
	hndl->init = init_solid_color;
	hndl->num_inputs = num_inputs_color;
	hndl->early_out = early_out_color;
	hndl->free = free_solid_color;
	hndl->copy = copy_solid_color;
	hndl->execute = do_solid_color;
}


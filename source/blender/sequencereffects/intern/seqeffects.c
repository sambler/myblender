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
 * - Shane Ambler 2011
 * - Shane Ambler 2011
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/sequencereffects/intern/seqeffects.c
 *  \ingroup seq
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "DNA_ID.h"
#include "DNA_sequence_types.h"

#include "BKE_sequencer.h"
#include "SEQ_effects.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

/* SeqConfig_* routines are each located in the files containing the effects code
 * -- only called within get_sequence_effect_impl() */

extern void SeqConfig_add(SeqEffect *hndl);				/* add.c			*/
extern void SeqConfig_adjustment(SeqEffect *hndl);		/* adjustment.c		*/
extern void SeqConfig_alphaOver(SeqEffect *hndl);		/* alpha.c			*/
extern void SeqConfig_alphaUnder(SeqEffect *hndl);		/* alpha.c			*/
extern void SeqConfig_cross(SeqEffect *hndl);			/* cross.c			*/
extern void SeqConfig_gamma(SeqEffect *hndl);			/* gamma.c			*/
extern void SeqConfig_glow(SeqEffect *hndl);			/* glow.c			*/
extern void SeqConfig_mul(SeqEffect *hndl);				/* mul.c			*/
extern void SeqConfig_multicam(SeqEffect *hndl);		/* multicam.c		*/
extern void SeqConfig_overdrop(SeqEffect *hndl);		/* drop.c			*/
extern void SeqConfig_plugin(SeqEffect *hndl);			/* plugin.c			*/
extern void SeqConfig_solid_colour(SeqEffect *hndl);	/* solid_colour.c	*/
extern void SeqConfig_speed(SeqEffect *hndl);			/* speed.c			*/
extern void SeqConfig_sub(SeqEffect *hndl);				/* sub.c			*/
extern void SeqConfig_transform(SeqEffect *hndl);		/* transform.c		*/
extern void SeqConfig_wipe(SeqEffect *hndl);			/* wipe.c			*/

static SeqEffect *fxList = NULL;

struct ImBuf * prepare_effect_imbufs(
	SeqRenderData context,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out;
	int x = context.rectx;
	int y = context.recty;

	if (!ibuf1 && !ibuf2 && !ibuf3) {
		/* hmmm, global float option ? */
		out = IMB_allocImBuf((short)x, (short)y, 32, IB_rect);
	} else if ((ibuf1 && ibuf1->rect_float) || 
				(ibuf2 && ibuf2->rect_float) ||
				(ibuf3 && ibuf3->rect_float)) {
		/* if any inputs are rectfloat, output is float too */

		out = IMB_allocImBuf((short)x, (short)y, 32, IB_rectfloat);
	} else {
		out = IMB_allocImBuf((short)x, (short)y, 32, IB_rect);
	}
	
	if (ibuf1 && !ibuf1->rect_float && out->rect_float) {
		IMB_float_from_rect_simple(ibuf1);
	}
	if (ibuf2 && !ibuf2->rect_float && out->rect_float) {
		IMB_float_from_rect_simple(ibuf2);
	}
	if (ibuf3 && !ibuf3->rect_float && out->rect_float) {
		IMB_float_from_rect_simple(ibuf3);
	}
	
	if (ibuf1 && !ibuf1->rect && !out->rect_float) {
		IMB_rect_from_float(ibuf1);
	}
	if (ibuf2 && !ibuf2->rect && !out->rect_float) {
		IMB_rect_from_float(ibuf2);
	}
	if (ibuf3 && !ibuf3->rect && !out->rect_float) {
		IMB_rect_from_float(ibuf3);
	}

	return out;
}



/* **********************************************************************
   sequence effect factory
   ********************************************************************** */

static void init_noop(struct Sequence *UNUSED(seq))
{

}

static void load_noop(struct Sequence *UNUSED(seq))
{

}

static void init_plugin_noop(struct Sequence *UNUSED(seq), const char *UNUSED(fname))
{

}

static void free_noop(struct Sequence *UNUSED(seq))
{

}

static int num_inputs_default(void)
{
	return 2;
}

static int early_out_noop(struct Sequence *UNUSED(seq),
			float UNUSED(facf0), float UNUSED(facf1))
{
	return 0;
}

static int early_out_fade(struct Sequence *UNUSED(seq),
			float facf0, float facf1)
{
	if (facf0 == 0.0f && facf1 == 0.0f) {
		return 1;
	} else if (facf0 == 1.0f && facf1 == 1.0f) {
		return 2;
	}
	return 0;
}

static int early_out_mul_input2(struct Sequence *UNUSED(seq),
				float facf0, float facf1)
{
	if (facf0 == 0.0f && facf1 == 0.0f) {
		return 1;
	}
	return 0;
}

static void store_icu_yrange_noop(struct Sequence * UNUSED(seq),
				short UNUSED(adrcode), float *UNUSED(ymin), float *UNUSED(ymax))
{
	/* defaults are fine */
}

static void get_default_fac_noop(struct Sequence *UNUSED(seq), float UNUSED(cfra),
				float * facf0, float * facf1)
{
	*facf0 = *facf1 = 1.0;
}

static void get_default_fac_fade(struct Sequence *seq, float cfra,
				float * facf0, float * facf1)
{
	*facf0 = (float)(cfra - seq->startdisp);
	*facf1 = (float)(*facf0 + 0.5f);
	*facf0 /= seq->len;
	*facf1 /= seq->len;
}

static SeqEffect get_sequence_effect_impl(int seq_type)
{
	SeqEffect rval;

	rval.init = init_noop;
	rval.init_plugin = init_plugin_noop;
	rval.num_inputs = num_inputs_default;
	rval.load = load_noop;
	rval.free = free_noop;
	rval.early_out = early_out_noop;
	rval.get_default_fac = get_default_fac_noop;
	rval.store_icu_yrange = store_icu_yrange_noop;
	rval.execute = NULL;
	rval.copy = NULL;

	switch (seq_type) {
	case SEQ_CROSS:
		rval.early_out = early_out_fade;
		rval.get_default_fac = get_default_fac_fade;
		SeqConfig_cross(&rval);
		break;
	case SEQ_GAMCROSS:
		rval.early_out = early_out_fade;
		rval.get_default_fac = get_default_fac_fade;
		SeqConfig_gamma(&rval);
		break;
	case SEQ_ADD:
		rval.early_out = early_out_mul_input2;
		SeqConfig_add(&rval);
		break;
	case SEQ_SUB:
		rval.early_out = early_out_mul_input2;
		SeqConfig_sub(&rval);
		break;
	case SEQ_MUL:
		rval.early_out = early_out_mul_input2;
		SeqConfig_mul(&rval);
		break;
	case SEQ_ALPHAOVER:
		SeqConfig_alphaOver(&rval);
		break;
	case SEQ_OVERDROP:
		SeqConfig_overdrop(&rval);
		break;
	case SEQ_ALPHAUNDER:
		SeqConfig_alphaUnder(&rval);
		break;
	case SEQ_WIPE:
		rval.early_out = early_out_fade;
		rval.get_default_fac = get_default_fac_fade;
		SeqConfig_wipe(&rval);
		break;
	case SEQ_GLOW:
		SeqConfig_glow(&rval);
		break;
	case SEQ_TRANSFORM:
		SeqConfig_transform(&rval);
		break;
	case SEQ_SPEED:
		SeqConfig_speed(&rval);
		break;
	case SEQ_COLOR:
		SeqConfig_solid_colour(&rval);
		break;
	case SEQ_PLUGIN:
		rval.get_default_fac = get_default_fac_fade;
		SeqConfig_plugin(&rval);
		break;
	case SEQ_MULTICAM:
		SeqConfig_multicam(&rval);
		break;
	case SEQ_ADJUSTMENT:
		SeqConfig_adjustment(&rval);
		break;
	}

	return rval;
}


SeqEffect get_sequence_effect(Sequence * seq)
{
	SeqEffect rval= {NULL};

	if (seq->type & SEQ_EFFECT) {
		rval = get_sequence_effect_impl(seq->type);
		if ((seq->flag & SEQ_EFFECT_NOT_LOADED) != 0) {
			rval.load(seq);
			seq->flag &= ~SEQ_EFFECT_NOT_LOADED;
		}
	}

	return rval;
}

SeqEffect get_sequence_blend(Sequence * seq)
{
	SeqEffect rval= {NULL};

	if (seq->blend_mode != 0) {
		rval = get_sequence_effect_impl(seq->blend_mode);
		if ((seq->flag & SEQ_EFFECT_NOT_LOADED) != 0) {
			rval.load(seq);
			seq->flag &= ~SEQ_EFFECT_NOT_LOADED;
		}
	}

	return rval;
}

int get_sequence_effect_num_inputs(int seq_type)
{
	SeqEffect rval = get_sequence_effect_impl(seq_type);

	int cnt = rval.num_inputs();
	if (rval.execute) {
		return cnt;
	}
	return 0;
}

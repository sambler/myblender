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
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/sequencereffects/intern/speed.c
 *  \ingroup seq
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"
#include "DNA_anim_types.h"

#include "BKE_sequencer.h"
#include "SEQ_effects.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

#include "RNA_access.h"

#include "seq_intern.h"

/* **********************************************************************
   SPEED
   ********************************************************************** */
static void init_speed_effect(Sequence *seq)
{
	SpeedControlVars * v;

	if(seq->effectdata) MEM_freeN(seq->effectdata);
	seq->effectdata = MEM_callocN(sizeof(struct SpeedControlVars), "speedcontrolvars");

	v = (SpeedControlVars *)seq->effectdata;
	v->globalSpeed = 1.0;
	v->frameMap = NULL;
	v->flags |= SEQ_SPEED_INTEGRATE; /* should be default behavior */
	v->length = 0;
}

static void load_speed_effect(Sequence * seq)
{
	SpeedControlVars * v = (SpeedControlVars *)seq->effectdata;

	v->frameMap = NULL;
	v->length = 0;
}

static int num_inputs_speed(void)
{
	return 1;
}

static void free_speed_effect(Sequence *seq)
{
	SpeedControlVars * v = (SpeedControlVars *)seq->effectdata;
	if(v->frameMap) MEM_freeN(v->frameMap);
	if(seq->effectdata) MEM_freeN(seq->effectdata);
	seq->effectdata = NULL;
}

static void copy_speed_effect(Sequence *dst, Sequence *src)
{
	SpeedControlVars * v;
	dst->effectdata = MEM_dupallocN(src->effectdata);
	v = (SpeedControlVars *)dst->effectdata;
	v->frameMap = NULL;
	v->length = 0;
}

static int early_out_speed(struct Sequence *UNUSED(seq), float UNUSED(facf0), float UNUSED(facf1))
{
	return 1;
}

static void store_icu_yrange_speed(struct Sequence * seq,
				short UNUSED(adrcode), float * ymin, float * ymax)
{
	SpeedControlVars * v = (SpeedControlVars *)seq->effectdata;

	/* if not already done, load / initialize data */
	get_sequence_effect(seq);

	if ((v->flags & SEQ_SPEED_INTEGRATE) != 0) {
		*ymin = -100.0;
		*ymax = 100.0;
	} else {
		if (v->flags & SEQ_SPEED_COMPRESS_IPO_Y) {
			*ymin = 0.0;
			*ymax = 1.0;
		} else {
			*ymin = 0.0;
			*ymax = seq->len;
		}
	}
}
void sequence_effect_speed_rebuild_map(Scene *scene, Sequence * seq, int force)
{
	int cfra;
	float fallback_fac = 1.0f;
	SpeedControlVars * v = (SpeedControlVars *)seq->effectdata;
	FCurve *fcu= NULL;
	int flags = v->flags;

	/* if not already done, load / initialize data */
	get_sequence_effect(seq);

	if ( (force == FALSE) && (seq->len == v->length) && (v->frameMap != NULL) ) {
		return;
	}
	if ( (seq->seq1 == NULL) || (seq->len < 1) ) {
		/* make coverity happy and check for (CID 598) input strip ... */
		return;
	}

	/* XXX - new in 2.5x. should we use the animation system this way?
	* The fcurve is needed because many frames need evaluating at once - campbell */
	fcu= id_data_find_fcurve(&scene->id, seq, &RNA_Sequence, "speed_factor", 0);


	if (!v->frameMap || v->length != seq->len) {
		if (v->frameMap) MEM_freeN(v->frameMap);

		v->length = seq->len;

		v->frameMap = MEM_callocN(sizeof(float) * v->length,
					"speedcontrol frameMap");
	}

	fallback_fac = 1.0;

	if (seq->flag & SEQ_USE_EFFECT_DEFAULT_FADE) {
		if (seq->seq1->enddisp != seq->seq1->start
		&& seq->seq1->len != 0) {
			fallback_fac = (float) seq->seq1->len /
				(float) (seq->seq1->enddisp - seq->seq1->start);
			flags = SEQ_SPEED_INTEGRATE;
			fcu = NULL;
		}
	} else {
		/* if there is no fcurve, use value as simple multiplier */
		if (!fcu) {
			fallback_fac = seq->speed_fader; /* same as speed_factor in rna*/
		}
	}

	if (flags & SEQ_SPEED_INTEGRATE) {
		float cursor = 0;
		float facf;

		v->frameMap[0] = 0;
		v->lastValidFrame = 0;

		for (cfra = 1; cfra < v->length; cfra++) {
			if(fcu) {
				facf = evaluate_fcurve(fcu, seq->startdisp + cfra);
			} else {
				facf = fallback_fac;
			}
			facf *= v->globalSpeed;

			cursor += facf;

			if (cursor >= seq->seq1->len) {
				v->frameMap[cfra] = seq->seq1->len - 1;
			} else {
				v->frameMap[cfra] = cursor;
				v->lastValidFrame = cfra;
			}
		}
	} else {
		float facf;

		v->lastValidFrame = 0;
		for (cfra = 0; cfra < v->length; cfra++) {

			if(fcu) {
				facf = evaluate_fcurve(fcu, seq->startdisp + cfra);
			} else {
				facf = fallback_fac;
			}

			if (flags & SEQ_SPEED_COMPRESS_IPO_Y) {
				facf *= seq->seq1->len;
			}
			facf *= v->globalSpeed;

			if (facf >= seq->seq1->len) {
				facf = seq->seq1->len - 1;
			} else {
				v->lastValidFrame = cfra;
			}
			v->frameMap[cfra] = facf;
		}
	}
}

/* setup */
void SeqConfigHandle_speed(struct SeqEffectHandle *hndl)
{
    hndl->init = init_speed_effect;
	hndl->num_inputs = num_inputs_speed;
	hndl->load = load_speed_effect;
	hndl->free = free_speed_effect;
	hndl->copy = copy_speed_effect;
	hndl->execute = do_cross_effect;
	hndl->early_out = early_out_speed;
    hndl->store_icu_yrange = store_icu_yrange_speed;
}


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

/** \file blender/sequencereffects/intern/cross.c
 *  \ingroup seq
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"

#include "BKE_sequencer.h"
#include "SEQ_effects.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

#include "seq_intern.h"

/* **********************************************************************
   CROSS
   ********************************************************************** */

static void do_cross_effect_byte(float facf0, float facf1, int x, int y,
			char *rect1, char *rect2, char *out)
{
	int fac1, fac2, fac3, fac4;
	int xo;
	char *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac2= (int)(256.0f*facf0);
	fac1= 256-fac2;
	fac4= (int)(256.0f*facf1);
	fac3= 256-fac4;

	while(y--) {

		x= xo;
		while(x--) {

			rt[0]= (fac1*rt1[0] + fac2*rt2[0])>>8;
			rt[1]= (fac1*rt1[1] + fac2*rt2[1])>>8;
			rt[2]= (fac1*rt1[2] + fac2*rt2[2])>>8;
			rt[3]= (fac1*rt1[3] + fac2*rt2[3])>>8;

			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			rt[0]= (fac3*rt1[0] + fac4*rt2[0])>>8;
			rt[1]= (fac3*rt1[1] + fac4*rt2[1])>>8;
			rt[2]= (fac3*rt1[2] + fac4*rt2[2])>>8;
			rt[3]= (fac3*rt1[3] + fac4*rt2[3])>>8;

			rt1+= 4; rt2+= 4; rt+= 4;
		}

	}
}

static void do_cross_effect_float(float facf0, float facf1, int x, int y,
			float *rect1, float *rect2, float *out)
{
	float fac1, fac2, fac3, fac4;
	int xo;
	float *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac2= facf0;
	fac1= 1.0f - fac2;
	fac4= facf1;
	fac3= 1.0f - fac4;

	while(y--) {

		x= xo;
		while(x--) {

			rt[0]= fac1*rt1[0] + fac2*rt2[0];
			rt[1]= fac1*rt1[1] + fac2*rt2[1];
			rt[2]= fac1*rt1[2] + fac2*rt2[2];
			rt[3]= fac1*rt1[3] + fac2*rt2[3];

			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			rt[0]= fac3*rt1[0] + fac4*rt2[0];
			rt[1]= fac3*rt1[1] + fac4*rt2[1];
			rt[2]= fac3*rt1[2] + fac4*rt2[2];
			rt[3]= fac3*rt1[3] + fac4*rt2[3];

			rt1+= 4; rt2+= 4; rt+= 4;
		}

	}
}

/* carefull: also used by speed effect! */

struct ImBuf* do_cross_effect(
	SeqRenderData context, Sequence *UNUSED(seq), float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(
		context, ibuf1, ibuf2, ibuf3);

	if (out->rect_float) {
		do_cross_effect_float(
			facf0, facf1, context.rectx, context.recty,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
	} else {
		do_cross_effect_byte(
			facf0, facf1, context.rectx, context.recty,
			(char*) ibuf1->rect, (char*) ibuf2->rect,
			(char*) out->rect);
	}
	return out;
}

/* setup */
void SeqConfigHandle_cross(struct SeqEffectHandle *hndl)
{
	hndl->execute = do_cross_effect;
}


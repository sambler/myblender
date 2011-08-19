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

/** \file blender/blenkernel/intern/seqeffects/seq_add.c
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
   ADD
   ********************************************************************** */

static void do_add_effect_byte(float facf0, float facf1, int x, int y,
				   unsigned char *rect1, unsigned char *rect2,
				   unsigned char *out)
{
	int col, xo, fac1, fac3;
	char *rt1, *rt2, *rt;

	xo= x;
	rt1= (char *)rect1;
	rt2= (char *)rect2;
	rt= (char *)out;

	fac1= (int)(256.0f*facf0);
	fac3= (int)(256.0f*facf1);

	while(y--) {

		x= xo;
		while(x--) {

			col= rt1[0]+ ((fac1*rt2[0])>>8);
			if(col>255) rt[0]= 255; else rt[0]= col;
			col= rt1[1]+ ((fac1*rt2[1])>>8);
			if(col>255) rt[1]= 255; else rt[1]= col;
			col= rt1[2]+ ((fac1*rt2[2])>>8);
			if(col>255) rt[2]= 255; else rt[2]= col;
			col= rt1[3]+ ((fac1*rt2[3])>>8);
			if(col>255) rt[3]= 255; else rt[3]= col;

			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			col= rt1[0]+ ((fac3*rt2[0])>>8);
			if(col>255) rt[0]= 255; else rt[0]= col;
			col= rt1[1]+ ((fac3*rt2[1])>>8);
			if(col>255) rt[1]= 255; else rt[1]= col;
			col= rt1[2]+ ((fac3*rt2[2])>>8);
			if(col>255) rt[2]= 255; else rt[2]= col;
			col= rt1[3]+ ((fac3*rt2[3])>>8);
			if(col>255) rt[3]= 255; else rt[3]= col;

			rt1+= 4; rt2+= 4; rt+= 4;
		}
	}
}

static void do_add_effect_float(float facf0, float facf1, int x, int y,
				float *rect1, float *rect2,
				float *out)
{
	int xo;
	float fac1, fac3;
	float *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac1= facf0;
	fac3= facf1;

	while(y--) {

		x= xo * 4;
		while(x--) {
			*rt = *rt1 + fac1 * (*rt2);

			rt1++; rt2++; rt++;
		}

		if(y==0) break;
		y--;

		x= xo * 4;
		while(x--) {
			*rt = *rt1 + fac3 * (*rt2);

			rt1++; rt2++; rt++;
		}
	}
}

static struct ImBuf * do_add_effect(SeqRenderData context,
				    Sequence *UNUSED(seq), float UNUSED(cfra),
				    float facf0, float facf1,
				    struct ImBuf *ibuf1, struct ImBuf *ibuf2,
				    struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	if (out->rect_float) {
		do_add_effect_float(
			facf0, facf1, context.rectx, context.recty,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
	} else {
		do_add_effect_byte(
			facf0, facf1, context.rectx, context.recty,
			(unsigned char*) ibuf1->rect, (unsigned char*) ibuf2->rect,
			(unsigned char*) out->rect);
	}
	return out;
}

/* setup */
void SeqConfigHandle_add(struct SeqEffectHandle *hndl)
{
	hndl->execute = do_add_effect;
}


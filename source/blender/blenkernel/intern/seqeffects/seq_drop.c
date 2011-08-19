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

/** \file blender/blenkernel/intern/seqeffects/seq_drop.c
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
   DROP
   ********************************************************************** */

/* Must be > 0 or add precopy, etc to the function */
#define XOFF	8
#define YOFF	8

static void do_drop_effect_byte(float facf0, float facf1, int x, int y,
				char *rect2i, char *rect1i, char *outi)
{
	int height, width, temp, fac, fac1, fac2;
	char *rt1, *rt2, *out;
	int field= 1;

	width= x;
	height= y;

	fac1= (int)(70.0f*facf0);
	fac2= (int)(70.0f*facf1);

	rt2= (char*) (rect2i + YOFF*width);
	rt1= (char*) rect1i;
	out= (char*) outi;
	for (y=0; y<height-YOFF; y++) {
		if(field) fac= fac1;
		else fac= fac2;
		field= !field;

		memcpy(out, rt1, sizeof(int)*XOFF);
		rt1+= XOFF*4;
		out+= XOFF*4;

		for (x=XOFF; x<width; x++) {
			temp= ((fac*rt2[3])>>8);

			*(out++)= MAX2(0, *rt1 - temp); rt1++;
			*(out++)= MAX2(0, *rt1 - temp); rt1++;
			*(out++)= MAX2(0, *rt1 - temp); rt1++;
			*(out++)= MAX2(0, *rt1 - temp); rt1++;
			rt2+=4;
		}
		rt2+=XOFF*4;
	}
	memcpy(out, rt1, sizeof(int)*YOFF*width);
}

static void do_drop_effect_float(float facf0, float facf1, int x, int y,
				float *rect2i, float *rect1i, float *outi)
{
	int height, width;
	float temp, fac, fac1, fac2;
	float *rt1, *rt2, *out;
	int field= 1;

	width= x;
	height= y;

	fac1= 70.0f*facf0;
	fac2= 70.0f*facf1;

	rt2=  (rect2i + YOFF*width);
	rt1=  rect1i;
	out=  outi;
	for (y=0; y<height-YOFF; y++) {
		if(field) fac= fac1;
		else fac= fac2;
		field= !field;

		memcpy(out, rt1, 4 * sizeof(float)*XOFF);
		rt1+= XOFF*4;
		out+= XOFF*4;

		for (x=XOFF; x<width; x++) {
			temp= fac * rt2[3];

			*(out++)= MAX2(0.0f, *rt1 - temp); rt1++;
			*(out++)= MAX2(0.0f, *rt1 - temp); rt1++;
			*(out++)= MAX2(0.0f, *rt1 - temp); rt1++;
			*(out++)= MAX2(0.0f, *rt1 - temp); rt1++;
			rt2+=4;
		}
		rt2+=XOFF*4;
	}
	memcpy(out, rt1, 4 * sizeof(float)*YOFF*width);
}

struct ImBuf * do_overdrop_effect(SeqRenderData context, Sequence *UNUSED(seq),
				float UNUSED(cfra), float facf0, float facf1,
				struct ImBuf * ibuf1, struct ImBuf * ibuf2, struct ImBuf * ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);
	int x = context.rectx;
	int y = context.recty;

	if (out->rect_float) {
		do_drop_effect_float(
			facf0, facf1, x, y,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
		do_alphaover_effect_float(
			facf0, facf1, x, y,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
	} else {
		do_drop_effect_byte(
			facf0, facf1, x, y,
			(char*) ibuf1->rect,
			(char*) ibuf2->rect,
			(char*) out->rect);
		do_alphaover_effect_byte(
			facf0, facf1, x, y,
			(char*) ibuf1->rect, (char*) ibuf2->rect,
			(char*) out->rect);
	}

	return out;
}


/* setup */
void SeqConfigHandle_overdrop(struct SeqEffectHandle *hndl)
{
	hndl->execute = do_overdrop_effect;
}


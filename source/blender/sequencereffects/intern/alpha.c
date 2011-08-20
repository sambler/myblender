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

/** \file blender/sequencereffects/intern/alpha.c
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
   ALPHA OVER
   ********************************************************************** */

static void init_alpha_over_or_under(Sequence * seq)
{
	Sequence * seq1 = seq->seq1;
	Sequence * seq2 = seq->seq2;

	seq->seq2= seq1;
	seq->seq1= seq2;
}

void do_alphaover_effect_byte(float facf0, float facf1, int x, int y,
					char * rect1, char *rect2, char *out)
{
	int fac2, mfac, fac, fac4;
	int xo, tempc;
	char *rt1, *rt2, *rt;

	xo= x;
	rt1= (char *)rect1;
	rt2= (char *)rect2;
	rt= (char *)out;

	fac2= (int)(256.0f*facf0);
	fac4= (int)(256.0f*facf1);

	while(y--) {

		x= xo;
		while(x--) {

			/* rt = rt1 over rt2  (alpha from rt1) */

			fac= fac2;
			mfac= 256 - ( (fac2*rt1[3])>>8 );

			if(fac==0) *( (unsigned int *)rt) = *( (unsigned int *)rt2);
			else if(mfac==0) *( (unsigned int *)rt) = *( (unsigned int *)rt1);
			else {
				tempc= ( fac*rt1[0] + mfac*rt2[0])>>8;
				if(tempc>255) rt[0]= 255; else rt[0]= tempc;
				tempc= ( fac*rt1[1] + mfac*rt2[1])>>8;
				if(tempc>255) rt[1]= 255; else rt[1]= tempc;
				tempc= ( fac*rt1[2] + mfac*rt2[2])>>8;
				if(tempc>255) rt[2]= 255; else rt[2]= tempc;
				tempc= ( fac*rt1[3] + mfac*rt2[3])>>8;
				if(tempc>255) rt[3]= 255; else rt[3]= tempc;
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			fac= fac4;
			mfac= 256 - ( (fac4*rt1[3])>>8 );

			if(fac==0) *( (unsigned int *)rt) = *( (unsigned int *)rt2);
			else if(mfac==0) *( (unsigned int *)rt) = *( (unsigned int *)rt1);
			else {
				tempc= ( fac*rt1[0] + mfac*rt2[0])>>8;
				if(tempc>255) rt[0]= 255; else rt[0]= tempc;
				tempc= ( fac*rt1[1] + mfac*rt2[1])>>8;
				if(tempc>255) rt[1]= 255; else rt[1]= tempc;
				tempc= ( fac*rt1[2] + mfac*rt2[2])>>8;
				if(tempc>255) rt[2]= 255; else rt[2]= tempc;
				tempc= ( fac*rt1[3] + mfac*rt2[3])>>8;
				if(tempc>255) rt[3]= 255; else rt[3]= tempc;
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}
	}
}

void do_alphaover_effect_float(float facf0, float facf1, int x, int y,
					float * rect1, float *rect2, float *out)
{
	float fac2, mfac, fac, fac4;
	int xo;
	float *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac2= facf0;
	fac4= facf1;

	while(y--) {

		x= xo;
		while(x--) {

			/* rt = rt1 over rt2  (alpha from rt1) */

			fac= fac2;
			mfac= 1.0f - (fac2*rt1[3]) ;

			if(fac <= 0.0f) {
				memcpy(rt, rt2, 4 * sizeof(float));
			} else if(mfac <=0) {
				memcpy(rt, rt1, 4 * sizeof(float));
			} else {
				rt[0] = fac*rt1[0] + mfac*rt2[0];
				rt[1] = fac*rt1[1] + mfac*rt2[1];
				rt[2] = fac*rt1[2] + mfac*rt2[2];
				rt[3] = fac*rt1[3] + mfac*rt2[3];
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			fac= fac4;
			mfac= 1.0f - (fac4*rt1[3]);

			if(fac <= 0.0f) {
				memcpy(rt, rt2, 4 * sizeof(float));
			} else if(mfac <= 0.0f) {
				memcpy(rt, rt1, 4 * sizeof(float));
			} else {
				rt[0] = fac*rt1[0] + mfac*rt2[0];
				rt[1] = fac*rt1[1] + mfac*rt2[1];
				rt[2] = fac*rt1[2] + mfac*rt2[2];
				rt[3] = fac*rt1[3] + mfac*rt2[3];
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}
	}
}

static struct ImBuf * do_alphaover_effect(
	SeqRenderData context, Sequence *UNUSED(seq), float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	if (out->rect_float) {
		do_alphaover_effect_float(
			facf0, facf1, context.rectx, context.recty,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
	} else {
		do_alphaover_effect_byte(
			facf0, facf1, context.rectx, context.recty,
			(char*) ibuf1->rect, (char*) ibuf2->rect,
			(char*) out->rect);
	}
	return out;
}


/* **********************************************************************
   ALPHA UNDER
   ********************************************************************** */

void do_alphaunder_effect_byte(
	float facf0, float facf1, int x, int y,
	char *rect1, char *rect2, char *out)
{
	int fac2, mfac, fac, fac4;
	int xo;
	char *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac2= (int)(256.0f*facf0);
	fac4= (int)(256.0f*facf1);

	while(y--) {

		x= xo;
		while(x--) {

			/* rt = rt1 under rt2  (alpha from rt2) */

			/* this complex optimalisation is because the
			 * 'skybuf' can be crossed in
			 */
			if(rt2[3]==0 && fac2==256) *( (unsigned int *)rt) = *( (unsigned int *)rt1);
			else if(rt2[3]==255) *( (unsigned int *)rt) = *( (unsigned int *)rt2);
			else {
				mfac= rt2[3];
				fac= (fac2*(256-mfac))>>8;

				if(fac==0) *( (unsigned int *)rt) = *( (unsigned int *)rt2);
				else {
					rt[0]= ( fac*rt1[0] + mfac*rt2[0])>>8;
					rt[1]= ( fac*rt1[1] + mfac*rt2[1])>>8;
					rt[2]= ( fac*rt1[2] + mfac*rt2[2])>>8;
					rt[3]= ( fac*rt1[3] + mfac*rt2[3])>>8;
				}
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			if(rt2[3]==0 && fac4==256) *( (unsigned int *)rt) = *( (unsigned int *)rt1);
			else if(rt2[3]==255) *( (unsigned int *)rt) = *( (unsigned int *)rt2);
			else {
				mfac= rt2[3];
				fac= (fac4*(256-mfac))>>8;

				if(fac==0) *( (unsigned int *)rt) = *( (unsigned int *)rt2);
				else {
					rt[0]= ( fac*rt1[0] + mfac*rt2[0])>>8;
					rt[1]= ( fac*rt1[1] + mfac*rt2[1])>>8;
					rt[2]= ( fac*rt1[2] + mfac*rt2[2])>>8;
					rt[3]= ( fac*rt1[3] + mfac*rt2[3])>>8;
				}
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}
	}
}


void do_alphaunder_effect_float(float facf0, float facf1, int x, int y,
					float *rect1, float *rect2, float *out)
{
	float fac2, mfac, fac, fac4;
	int xo;
	float *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac2= facf0;
	fac4= facf1;

	while(y--) {

		x= xo;
		while(x--) {

			/* rt = rt1 under rt2  (alpha from rt2) */

			/* this complex optimalisation is because the
			 * 'skybuf' can be crossed in
			 */
			if( rt2[3]<=0 && fac2 >= 1.0f) {
				memcpy(rt, rt1, 4 * sizeof(float));
			} else if(rt2[3] >= 1.0f) {
				memcpy(rt, rt2, 4 * sizeof(float));
			} else {
				mfac = rt2[3];
				fac = fac2 * (1.0f - mfac);

				if(fac == 0) {
					memcpy(rt, rt2, 4 * sizeof(float));
				} else {
					rt[0]= fac*rt1[0] + mfac*rt2[0];
					rt[1]= fac*rt1[1] + mfac*rt2[1];
					rt[2]= fac*rt1[2] + mfac*rt2[2];
					rt[3]= fac*rt1[3] + mfac*rt2[3];
				}
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			if(rt2[3]<=0 && fac4 >= 1.0f) {
				memcpy(rt, rt1, 4 * sizeof(float));

			} else if(rt2[3]>=1.0f) {
				memcpy(rt, rt2, 4 * sizeof(float));
			} else {
				mfac= rt2[3];
				fac= fac4*(1.0f-mfac);

				if(fac == 0) {
					memcpy(rt, rt2, 4 * sizeof(float));
				} else {
					rt[0]= fac * rt1[0] + mfac * rt2[0];
					rt[1]= fac * rt1[1] + mfac * rt2[1];
					rt[2]= fac * rt1[2] + mfac * rt2[2];
					rt[3]= fac * rt1[3] + mfac * rt2[3];
				}
			}
			rt1+= 4; rt2+= 4; rt+= 4;
		}
	}
}

static struct ImBuf* do_alphaunder_effect(
	SeqRenderData context, Sequence *UNUSED(seq), float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(
		context, ibuf1, ibuf2, ibuf3);

	if (out->rect_float) {
		do_alphaunder_effect_float(
			facf0, facf1, context.rectx, context.recty,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
	} else {
		do_alphaunder_effect_byte(
			facf0, facf1, context.rectx, context.recty,
			(char*) ibuf1->rect, (char*) ibuf2->rect,
			(char*) out->rect);
	}
	return out;
}

/* setup */
void SeqConfigHandle_alphaOver(struct SeqEffectHandle *hndl)
{
	hndl->init = init_alpha_over_or_under;
	hndl->execute = do_alphaover_effect;
}

void SeqConfigHandle_alphaUnder(struct SeqEffectHandle *hndl)
{
	hndl->init = init_alpha_over_or_under;
	hndl->execute = do_alphaunder_effect;
}


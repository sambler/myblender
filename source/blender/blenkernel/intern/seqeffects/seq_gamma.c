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

/** \file blender/blenkernel/intern/seqeffects/seq_gamma.c
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
   GAMMA CROSS
   ********************************************************************** */

/* copied code from initrender.c */
static unsigned short gamtab[65536];
static unsigned short igamtab1[256];
static int gamma_tabs_init = FALSE;

#define RE_GAMMA_TABLE_SIZE 400

static float gamma_range_table[RE_GAMMA_TABLE_SIZE + 1];
static float gamfactor_table[RE_GAMMA_TABLE_SIZE];
static float inv_gamma_range_table[RE_GAMMA_TABLE_SIZE + 1];
static float inv_gamfactor_table[RE_GAMMA_TABLE_SIZE];
static float color_domain_table[RE_GAMMA_TABLE_SIZE + 1];
static float color_step;
static float inv_color_step;
static float valid_gamma;
static float valid_inv_gamma;

static void makeGammaTables(float gamma)
{
	/* we need two tables: one forward, one backward */
	int i;

	valid_gamma        = gamma;
	valid_inv_gamma    = 1.0f / gamma;
	color_step        = 1.0f / RE_GAMMA_TABLE_SIZE;
	inv_color_step    = (float) RE_GAMMA_TABLE_SIZE;

	/* We could squeeze out the two range tables to gain some memory.        */
	for (i = 0; i < RE_GAMMA_TABLE_SIZE; i++) {
		color_domain_table[i]   = i * color_step;
		gamma_range_table[i]     = pow(color_domain_table[i],
										valid_gamma);
		inv_gamma_range_table[i] = pow(color_domain_table[i],
										valid_inv_gamma);
	}

	/* The end of the table should match 1.0 carefully. In order to avoid    */
	/* rounding errors, we just set this explicitly. The last segment may    */
	/* have a different length than the other segments, but our              */
	/* interpolation is insensitive to that.                                 */
	color_domain_table[RE_GAMMA_TABLE_SIZE]   = 1.0;
	gamma_range_table[RE_GAMMA_TABLE_SIZE]     = 1.0;
	inv_gamma_range_table[RE_GAMMA_TABLE_SIZE] = 1.0;

	/* To speed up calculations, we make these calc factor tables. They are  */
	/* multiplication factors used in scaling the interpolation.             */
	for (i = 0; i < RE_GAMMA_TABLE_SIZE; i++ ) {
		gamfactor_table[i] = inv_color_step
			* (gamma_range_table[i + 1] - gamma_range_table[i]) ;
		inv_gamfactor_table[i] = inv_color_step
			* (inv_gamma_range_table[i + 1] - inv_gamma_range_table[i]) ;
	}

} /* end of void makeGammaTables(float gamma) */


static float gammaCorrect(float c)
{
	int i;
	float res = 0.0;

	i = floor(c * inv_color_step);
	/* Clip to range [0,1]: outside, just do the complete calculation.       */
	/* We may have some performance problems here. Stretching up the LUT     */
	/* may help solve that, by exchanging LUT size for the interpolation.    */
	/* Negative colors are explicitly handled.                              */
	if (i < 0) res = -pow(abs(c), valid_gamma);
	else if (i >= RE_GAMMA_TABLE_SIZE ) res = pow(c, valid_gamma);
	else res = gamma_range_table[i] +
			   ( (c - color_domain_table[i]) * gamfactor_table[i]);

	return res;
} /* end of float gammaCorrect(float col) */

/* ------------------------------------------------------------------------- */

static float invGammaCorrect(float col)
{
	int i;
	float res = 0.0;

	i = floor(col*inv_color_step);
	/* Negative colors are explicitly handled.                              */
	if (i < 0) res = -pow(abs(col), valid_inv_gamma);
	else if (i >= RE_GAMMA_TABLE_SIZE) res = pow(col, valid_inv_gamma);
	else res = inv_gamma_range_table[i] +
			   ( (col - color_domain_table[i]) * inv_gamfactor_table[i]);

	return res;
} /* end of float invGammaCorrect(float col) */


static void gamtabs(float gamma)
{
	float val, igamma= 1.0f/gamma;
	int a;

	/* gamtab: in short, out short */
	for(a=0; a<65536; a++) {
		val= a;
		val/= 65535.0f;

		if(gamma==2.0f) val= sqrt(val);
		else if(gamma!=1.0f) val= pow(val, igamma);

		gamtab[a]= (65535.99f*val);
	}
	/* inverse gamtab1 : in byte, out short */
	for(a=1; a<=256; a++) {
		if(gamma==2.0f) igamtab1[a-1]= a*a-1;
		else if(gamma==1.0f) igamtab1[a-1]= 256*a-1;
		else {
			val= a/256.0f;
			igamtab1[a-1]= (65535.0*pow(val, gamma)) -1 ;
		}
	}

}

static void build_gammatabs(void)
{
	if (gamma_tabs_init == FALSE) {
		gamtabs(2.0f);
		makeGammaTables(2.0f);
		gamma_tabs_init = TRUE;
	}
}

static void init_gammacross(Sequence * UNUSED(seq))
{
}

static void load_gammacross(Sequence * UNUSED(seq))
{
}

static void free_gammacross(Sequence * UNUSED(seq))
{
}

static void do_gammacross_effect_byte(float facf0, float UNUSED(facf1),
					  int x, int y,
					  unsigned char *rect1,
					  unsigned char *rect2,
					  unsigned char *out)
{
	int fac1, fac2, col;
	int xo;
	unsigned char *rt1, *rt2, *rt;

	xo= x;
	rt1= (unsigned char *)rect1;
	rt2= (unsigned char *)rect2;
	rt= (unsigned char *)out;

	fac2= (int)(256.0f*facf0);
	fac1= 256-fac2;

	while(y--) {

		x= xo;
		while(x--) {

			col= (fac1*igamtab1[rt1[0]] + fac2*igamtab1[rt2[0]])>>8;
			if(col>65535) rt[0]= 255; else rt[0]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];
			col=(fac1*igamtab1[rt1[1]] + fac2*igamtab1[rt2[1]])>>8;
			if(col>65535) rt[1]= 255; else rt[1]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];
			col= (fac1*igamtab1[rt1[2]] + fac2*igamtab1[rt2[2]])>>8;
			if(col>65535) rt[2]= 255; else rt[2]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];
			col= (fac1*igamtab1[rt1[3]] + fac2*igamtab1[rt2[3]])>>8;
			if(col>65535) rt[3]= 255; else rt[3]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];

			rt1+= 4; rt2+= 4; rt+= 4;
		}

		if(y==0) break;
		y--;

		x= xo;
		while(x--) {

			col= (fac1*igamtab1[rt1[0]] + fac2*igamtab1[rt2[0]])>>8;
			if(col>65535) rt[0]= 255; else rt[0]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];
			col= (fac1*igamtab1[rt1[1]] + fac2*igamtab1[rt2[1]])>>8;
			if(col>65535) rt[1]= 255; else rt[1]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];
			col= (fac1*igamtab1[rt1[2]] + fac2*igamtab1[rt2[2]])>>8;
			if(col>65535) rt[2]= 255; else rt[2]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];
			col= (fac1*igamtab1[rt1[3]] + fac2*igamtab1[rt2[3]])>>8;
			if(col>65535) rt[3]= 255; else rt[3]= ( (char *)(gamtab+col))[MOST_SIG_BYTE];

			rt1+= 4; rt2+= 4; rt+= 4;
		}
	}

}

static void do_gammacross_effect_float(float facf0, float UNUSED(facf1),
					   int x, int y,
					   float *rect1, float *rect2,
					   float *out)
{
	float fac1, fac2;
	int xo;
	float *rt1, *rt2, *rt;

	xo= x;
	rt1= rect1;
	rt2= rect2;
	rt= out;

	fac2= facf0;
	fac1= 1.0f - fac2;

	while(y--) {

		x= xo * 4;
		while(x--) {

			*rt= gammaCorrect(
				fac1 * invGammaCorrect(*rt1)
				+ fac2 * invGammaCorrect(*rt2));
			rt1++; rt2++; rt++;
		}

		if(y==0) break;
		y--;

		x= xo * 4;
		while(x--) {

			*rt= gammaCorrect(
				fac1*invGammaCorrect(*rt1)
				+ fac2*invGammaCorrect(*rt2));

			rt1++; rt2++; rt++;
		}
	}
}

static struct ImBuf * do_gammacross_effect(
	SeqRenderData context,
	Sequence *UNUSED(seq), float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	build_gammatabs();

	if (out->rect_float) {
		do_gammacross_effect_float(
			facf0, facf1, context.rectx, context.recty,
			ibuf1->rect_float, ibuf2->rect_float,
			out->rect_float);
	} else {
		do_gammacross_effect_byte(
			facf0, facf1, context.rectx, context.recty,
			(unsigned char*) ibuf1->rect, (unsigned char*) ibuf2->rect,
			(unsigned char*) out->rect);
	}
	return out;
}

/* setup */
void SeqConfigHandle_gamma(struct SeqEffectHandle *hndl)
{
	hndl->init = init_gammacross;
	hndl->load = load_gammacross;
	hndl->free = free_gammacross;
	hndl->execute = do_gammacross_effect;
}


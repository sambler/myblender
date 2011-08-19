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

/** \file blender/blenkernel/intern/seqeffects/seq_wipe.c
 *  \ingroup bke
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_math.h" /* windows needs for M_PI */
#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"

#include "BKE_sequencer.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

#include "seq_intern.h"

/* **********************************************************************
   WIPE
   ********************************************************************** */

typedef struct WipeZone {
	float angle;
	int flip;
	int xo, yo;
	int width;
	float invwidth;
	float pythangle;
} WipeZone;

static void precalc_wipe_zone(WipeZone *wipezone, WipeVars *wipe, int xo, int yo)
{
	wipezone->flip = (wipe->angle < 0);
	wipezone->angle = pow(fabsf(wipe->angle)/45.0f, log(xo)/M_LN2);
	wipezone->xo = xo;
	wipezone->yo = yo;
	wipezone->width = (int)(wipe->edgeWidth*((xo+yo)/2.0f));
	wipezone->pythangle = 1.0f/sqrtf(wipe->angle*wipe->angle + 1.0f);

	if(wipe->wipetype == DO_SINGLE_WIPE)
		wipezone->invwidth = 1.0f/wipezone->width;
	else
		wipezone->invwidth = 1.0f/(0.5f*wipezone->width);
}

// This function calculates the blur band for the wipe effects
static float in_band(WipeZone *wipezone,float width,float dist,float perc,int side,int dir)
{
	float t1,t2,alpha;

	if(width == 0)
		return (float)side;

	if(width < dist)
		return side;

	t1 = dist * wipezone->invwidth;  //percentange of width that is
	t2 = wipezone->invwidth;  //amount of alpha per % point

	if(side == 1)
		alpha = (t1*t2*100) + (1-perc); // add point's alpha contrib to current position in wipe
	else
		alpha = (1-perc) - (t1*t2*100);

	if(dir == 0)
		alpha = 1-alpha;

	return alpha;
}

static float check_zone(WipeZone *wipezone, int x, int y,
	Sequence *seq, float facf0)
{
	float posx, posy,hyp,hyp2,angle,hwidth,b1,b2,b3,pointdist;
/*some future stuff
float hyp3,hyp4,b4,b5
*/
	float temp1,temp2,temp3,temp4; //some placeholder variables
	int xo = wipezone->xo;
	int yo = wipezone->yo;
	float halfx = xo*0.5f;
	float halfy = yo*0.5f;
	float widthf,output=0;
	WipeVars *wipe = (WipeVars *)seq->effectdata;
	int width;

	if(wipezone->flip) x = xo - x;
	angle = wipezone->angle;

	posy = facf0 * yo;

	if(wipe->forward){
		posx = facf0 * xo;
		posy = facf0 * yo;
	} else{
		posx = xo - facf0 * xo;
		posy = yo - facf0 * yo;
	}

	switch (wipe->wipetype) {
		case DO_SINGLE_WIPE:
			width = wipezone->width;
			hwidth = width*0.5f;

			if(angle == 0.0f) {
				b1 = posy;
				b2 = y;
				hyp = fabs(y - posy);
			}
			else {
				b1 = posy - (-angle)*posx;
				b2 = y - (-angle)*x;
				hyp = fabsf(angle*x+y+(-posy-angle*posx))*wipezone->pythangle;
			}

			if(angle < 0) {
				temp1 = b1;
				b1 = b2;
				b2 = temp1;
			}

			if(wipe->forward) {
				if(b1 < b2)
					output = in_band(wipezone,width,hyp,facf0,1,1);
				else
					output = in_band(wipezone,width,hyp,facf0,0,1);
			}
			else {
				if(b1 < b2)
					output = in_band(wipezone,width,hyp,facf0,0,1);
				else
					output = in_band(wipezone,width,hyp,facf0,1,1);
			}
		break;

		case DO_DOUBLE_WIPE:
			if(!wipe->forward)
				facf0 = 1.0f-facf0;   // Go the other direction

			width = wipezone->width;  // calculate the blur width
			hwidth = width*0.5f;
			if (angle == 0) {
				b1 = posy*0.5f;
				b3 = yo-posy*0.5f;
				b2 = y;

				hyp = abs(y - posy*0.5f);
				hyp2 = abs(y - (yo-posy*0.5f));
			}
			else {
				b1 = posy*0.5f - (-angle)*posx*0.5f;
				b3 = (yo-posy*0.5f) - (-angle)*(xo-posx*0.5f);
				b2 = y - (-angle)*x;

				hyp = abs(angle*x+y+(-posy*0.5f-angle*posx*0.5f))*wipezone->pythangle;
				hyp2 = abs(angle*x+y+(-(yo-posy*0.5f)-angle*(xo-posx*0.5f)))*wipezone->pythangle;
			}

			temp1 = xo*(1-facf0*0.5f)-xo*facf0*0.5f;
			temp2 = yo*(1-facf0*0.5f)-yo*facf0*0.5f;
			pointdist = sqrt(temp1*temp1 + temp2*temp2);

			if(b2 < b1 && b2 < b3 ){
				if(hwidth < pointdist)
					output = in_band(wipezone,hwidth,hyp,facf0,0,1);
			} else if(b2 > b1 && b2 > b3 ){
				if(hwidth < pointdist)
					output = in_band(wipezone,hwidth,hyp2,facf0,0,1);
			} else {
				if(  hyp < hwidth && hyp2 > hwidth )
					output = in_band(wipezone,hwidth,hyp,facf0,1,1);
				else if( hyp > hwidth && hyp2 < hwidth )
					  output = in_band(wipezone,hwidth,hyp2,facf0,1,1);
				else
					  output = in_band(wipezone,hwidth,hyp2,facf0,1,1) * in_band(wipezone,hwidth,hyp,facf0,1,1);
			}
			if(!wipe->forward)output = 1-output;
		break;
		case DO_CLOCK_WIPE:
			  /*
				  temp1: angle of effect center in rads
				  temp2: angle of line through (halfx,halfy) and (x,y) in rads
				  temp3: angle of low side of blur
				  temp4: angle of high side of blur
			  */
			output = 1.0f - facf0;
			widthf = wipe->edgeWidth*2.0f*(float)M_PI;
			temp1 = 2.0f * (float)M_PI * facf0;

			if(wipe->forward){
				temp1 = 2.0f*(float)M_PI - temp1;
			}

			x = x - halfx;
			y = y - halfy;

			temp2 = asin(abs(y)/sqrt(x*x + y*y));
			if(x <= 0 && y >= 0) temp2 = (float)M_PI - temp2;
			else if(x<=0 && y <= 0) temp2 += (float)M_PI;
			else if(x >= 0 && y <= 0) temp2 = 2.0f*(float)M_PI - temp2;

			if(wipe->forward){
				temp3 = temp1-(widthf*0.5f)*facf0;
				temp4 = temp1+(widthf*0.5f)*(1-facf0);
			} else{
				temp3 = temp1-(widthf*0.5f)*(1-facf0);
				temp4 = temp1+(widthf*0.5f)*facf0;
			}
			if (temp3 < 0) temp3 = 0;
			if (temp4 > 2.0f*(float)M_PI) temp4 = 2.0f*(float)M_PI;


			if(temp2 < temp3) output = 0;
			else if (temp2 > temp4) output = 1;
			else output = (temp2-temp3)/(temp4-temp3);
			if(x == 0 && y == 0) output = 1;
			if(output != output) output = 1;
			if(wipe->forward) output = 1 - output;
		break;
	/* BOX WIPE IS NOT WORKING YET */
	/* case DO_CROSS_WIPE: */
	/* BOX WIPE IS NOT WORKING YET */
	/*
		case DO_BOX_WIPE:
			if(invert)facf0 = 1-facf0;

			width = (int)(wipe->edgeWidth*((xo+yo)/2.0));
			hwidth = (float)width/2.0;
			if (angle == 0)angle = 0.000001;
			b1 = posy/2 - (-angle)*posx/2;
			b3 = (yo-posy/2) - (-angle)*(xo-posx/2);
			b2 = y - (-angle)*x;

			hyp = abs(angle*x+y+(-posy/2-angle*posx/2))*wipezone->pythangle;
			hyp2 = abs(angle*x+y+(-(yo-posy/2)-angle*(xo-posx/2)))*wipezone->pythangle;

			temp1 = xo*(1-facf0/2)-xo*facf0/2;
			temp2 = yo*(1-facf0/2)-yo*facf0/2;
			pointdist = sqrt(temp1*temp1 + temp2*temp2);

			if(b2 < b1 && b2 < b3 ){
				if(hwidth < pointdist)
					output = in_band(wipezone,hwidth,hyp,facf0,0,1);
			} else if(b2 > b1 && b2 > b3 ){
				if(hwidth < pointdist)
					output = in_band(wipezone,hwidth,hyp2,facf0,0,1);
			} else {
				if( hyp < hwidth && hyp2 > hwidth )
					output = in_band(wipezone,hwidth,hyp,facf0,1,1);
				else if( hyp > hwidth && hyp2 < hwidth )
					 output = in_band(wipezone,hwidth,hyp2,facf0,1,1);
				else
					 output = in_band(wipezone,hwidth,hyp2,facf0,1,1) * in_band(wipezone,hwidth,hyp,facf0,1,1);
			}

			if(invert)facf0 = 1-facf0;
			angle = -1/angle;
			b1 = posy/2 - (-angle)*posx/2;
			b3 = (yo-posy/2) - (-angle)*(xo-posx/2);
			b2 = y - (-angle)*x;

			hyp = abs(angle*x+y+(-posy/2-angle*posx/2))*wipezone->pythangle;
			hyp2 = abs(angle*x+y+(-(yo-posy/2)-angle*(xo-posx/2)))*wipezone->pythangle;

			if(b2 < b1 && b2 < b3 ){
				if(hwidth < pointdist)
					output *= in_band(wipezone,hwidth,hyp,facf0,0,1);
			} else if(b2 > b1 && b2 > b3 ){
				if(hwidth < pointdist)
					output *= in_band(wipezone,hwidth,hyp2,facf0,0,1);
			} else {
				if( hyp < hwidth && hyp2 > hwidth )
					output *= in_band(wipezone,hwidth,hyp,facf0,1,1);
				else if( hyp > hwidth && hyp2 < hwidth )
					output *= in_band(wipezone,hwidth,hyp2,facf0,1,1);
				else
					output *= in_band(wipezone,hwidth,hyp2,facf0,1,1) * in_band(wipezone,hwidth,hyp,facf0,1,1);
			}

		break;
*/
		case DO_IRIS_WIPE:
			if(xo > yo) yo = xo;
			else xo = yo;

			if(!wipe->forward) facf0 = 1-facf0;

			width = wipezone->width;
			hwidth = width*0.5f;

			temp1 = (halfx-(halfx)*facf0);
			pointdist = sqrt(temp1*temp1 + temp1*temp1);

			temp2 = sqrt((halfx-x)*(halfx-x) + (halfy-y)*(halfy-y));
			if(temp2 > pointdist) output = in_band(wipezone,hwidth,fabs(temp2-pointdist),facf0,0,1);
			else output = in_band(wipezone,hwidth,fabs(temp2-pointdist),facf0,1,1);

			if(!wipe->forward) output = 1-output;

		break;
	}
	if (output < 0) output = 0;
	else if(output > 1) output = 1;
	return output;
}

static void init_wipe_effect(Sequence *seq)
{
	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = MEM_callocN(sizeof(struct WipeVars), "wipevars");
}

static int num_inputs_wipe(void)
{
	return 1;
}

static void free_wipe_effect(Sequence *seq)
{
	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = NULL;
}

static void copy_wipe_effect(Sequence *dst, Sequence *src)
{
	dst->effectdata = MEM_dupallocN(src->effectdata);
}

static void do_wipe_effect_byte(Sequence *seq, float facf0, float UNUSED(facf1),
				int x, int y,
				unsigned char *rect1,
				unsigned char *rect2, unsigned char *out)
{
	WipeZone wipezone;
	WipeVars *wipe = (WipeVars *)seq->effectdata;
	int xo, yo;
	char *rt1, *rt2, *rt;

	precalc_wipe_zone(&wipezone, wipe, x, y);

	rt1 = (char *)rect1;
	rt2 = (char *)rect2;
	rt = (char *)out;

	xo = x;
	yo = y;
	for(y=0;y<yo;y++) {
		for(x=0;x<xo;x++) {
			float check = check_zone(&wipezone,x,y,seq,facf0);
			if (check) {
				if (rt1) {
					rt[0] = (int)(rt1[0]*check)+ (int)(rt2[0]*(1-check));
					rt[1] = (int)(rt1[1]*check)+ (int)(rt2[1]*(1-check));
					rt[2] = (int)(rt1[2]*check)+ (int)(rt2[2]*(1-check));
					rt[3] = (int)(rt1[3]*check)+ (int)(rt2[3]*(1-check));
				} else {
					rt[0] = 0;
					rt[1] = 0;
					rt[2] = 0;
					rt[3] = 255;
				}
			} else {
				if (rt2) {
					rt[0] = rt2[0];
					rt[1] = rt2[1];
					rt[2] = rt2[2];
					rt[3] = rt2[3];
				} else {
					rt[0] = 0;
					rt[1] = 0;
					rt[2] = 0;
					rt[3] = 255;
				}
			}

			rt+=4;
			if(rt1 !=NULL){
				rt1+=4;
			}
			if(rt2 !=NULL){
				rt2+=4;
			}
		}
	}
}

static void do_wipe_effect_float(Sequence *seq, float facf0, float UNUSED(facf1),
				 int x, int y,
				 float *rect1,
				 float *rect2, float *out)
{
	WipeZone wipezone;
	WipeVars *wipe = (WipeVars *)seq->effectdata;
	int xo, yo;
	float *rt1, *rt2, *rt;

	precalc_wipe_zone(&wipezone, wipe, x, y);

	rt1 = rect1;
	rt2 = rect2;
	rt = out;

	xo = x;
	yo = y;
	for(y=0;y<yo;y++) {
		for(x=0;x<xo;x++) {
			float check = check_zone(&wipezone,x,y,seq,facf0);
			if (check) {
				if (rt1) {
					rt[0] = rt1[0]*check+ rt2[0]*(1-check);
					rt[1] = rt1[1]*check+ rt2[1]*(1-check);
					rt[2] = rt1[2]*check+ rt2[2]*(1-check);
					rt[3] = rt1[3]*check+ rt2[3]*(1-check);
				} else {
					rt[0] = 0;
					rt[1] = 0;
					rt[2] = 0;
					rt[3] = 1.0;
				}
			} else {
				if (rt2) {
					rt[0] = rt2[0];
					rt[1] = rt2[1];
					rt[2] = rt2[2];
					rt[3] = rt2[3];
				} else {
					rt[0] = 0;
					rt[1] = 0;
					rt[2] = 0;
					rt[3] = 1.0;
				}
			}

			rt+=4;
			if(rt1 !=NULL){
				rt1+=4;
			}
			if(rt2 !=NULL){
				rt2+=4;
			}
		}
	}
}

static struct ImBuf * do_wipe_effect(
	SeqRenderData context, Sequence *seq, float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	if (out->rect_float) {
		do_wipe_effect_float(seq,
				     facf0, facf1, context.rectx, context.recty,
				     ibuf1->rect_float, ibuf2->rect_float,
				     out->rect_float);
	} else {
		do_wipe_effect_byte(seq,
				    facf0, facf1, context.rectx, context.recty,
				    (unsigned char*) ibuf1->rect, (unsigned char*) ibuf2->rect,
				    (unsigned char*) out->rect);
	}

	return out;
}

/* setup */
void SeqConfigHandle_wipe(struct SeqEffectHandle *hndl)
{
	hndl->init = init_wipe_effect;
	hndl->num_inputs = num_inputs_wipe;
	hndl->free = free_wipe_effect;
	hndl->copy = copy_wipe_effect;
	hndl->execute = do_wipe_effect;
}


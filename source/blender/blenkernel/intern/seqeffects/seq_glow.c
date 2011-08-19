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

/** \file blender/blenkernel/intern/seqeffects/seq_glow.c
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

enum {
	GlowR=0,
	GlowG=1,
	GlowB=2,
	GlowA=3
};

/* **********************************************************************
   GLOW
   ********************************************************************** */

static void RVBlurBitmap2_byte ( unsigned char* map, int width,int height,
				float blur, int quality)
/*	MUUUCCH better than the previous blur. */
/*	We do the blurring in two passes which is a whole lot faster. */
/*	I changed the math arount to implement an actual Gaussian */
/*	distribution. */
/*	*/
/*	Watch out though, it tends to misbehaven with large blur values on */
/*	a small bitmap.  Avoid avoid avoid. */
/*=============================== */
{
	unsigned char*	temp=NULL,*swap;
	float	*filter=NULL;
	int		x,y,i,fx,fy;
	int		index, ix, halfWidth;
	float	fval, k, curColor[3], curColor2[3], weight=0;

	/*	If we're not really blurring, bail out */
	if (blur<=0)
		return;

	/*	Allocate memory for the tempmap and the blur filter matrix */
	temp= MEM_mallocN( (width*height*4), "blurbitmaptemp");
	if (!temp)
		return;

	/*	Allocate memory for the filter elements */
	halfWidth = ((quality+1)*blur);
	filter = (float *)MEM_mallocN(sizeof(float)*halfWidth*2, "blurbitmapfilter");
	if (!filter){
		MEM_freeN (temp);
		return;
	}

	/*	Apparently we're calculating a bell curve */
	/*	based on the standard deviation (or radius) */
	/*	This code is based on an example */
	/*	posted to comp.graphics.algorithms by */
	/*	Blancmange (bmange@airdmhor.gen.nz) */

	k = -1.0f/(2.0f*(float)M_PI*blur*blur);
	for (ix = 0;ix< halfWidth;ix++){
		weight = (float)exp(k*(ix*ix));
		filter[halfWidth - ix] = weight;
		filter[halfWidth + ix] = weight;
	}
	filter[0] = weight;

	/*	Normalize the array */
	fval=0;
	for (ix = 0;ix< halfWidth*2;ix++)
		fval+=filter[ix];

	for (ix = 0;ix< halfWidth*2;ix++)
		filter[ix]/=fval;

	/*	Blur the rows */
	for (y=0;y<height;y++){
		/*	Do the left & right strips */
		for (x=0;x<halfWidth;x++){
			index=(x+y*width)*4;
			fx=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			curColor2[0]=curColor2[1]=curColor2[2]=0;

			for (i=x-halfWidth;i<x+halfWidth;i++){
				if ((i>=0)&&(i<width)){
					curColor[0]+=map[(i+y*width)*4+GlowR]*filter[fx];
					curColor[1]+=map[(i+y*width)*4+GlowG]*filter[fx];
					curColor[2]+=map[(i+y*width)*4+GlowB]*filter[fx];

					curColor2[0]+=map[(width-1-i+y*width)*4+GlowR] *
						filter[fx];
					curColor2[1]+=map[(width-1-i+y*width)*4+GlowG] *
						filter[fx];
					curColor2[2]+=map[(width-1-i+y*width)*4+GlowB] *
						filter[fx];
				}
				fx++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];

			temp[((width-1-x+y*width)*4)+GlowR]=curColor2[0];
			temp[((width-1-x+y*width)*4)+GlowG]=curColor2[1];
			temp[((width-1-x+y*width)*4)+GlowB]=curColor2[2];

		}
		/*	Do the main body */
		for (x=halfWidth;x<width-halfWidth;x++){
			index=(x+y*width)*4;
			fx=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			for (i=x-halfWidth;i<x+halfWidth;i++){
				curColor[0]+=map[(i+y*width)*4+GlowR]*filter[fx];
				curColor[1]+=map[(i+y*width)*4+GlowG]*filter[fx];
				curColor[2]+=map[(i+y*width)*4+GlowB]*filter[fx];
				fx++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];
		}
	}

	/*	Swap buffers */
	swap=temp;temp=map;map=swap;


	/*	Blur the columns */
	for (x=0;x<width;x++){
		/*	Do the top & bottom strips */
		for (y=0;y<halfWidth;y++){
			index=(x+y*width)*4;
			fy=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			curColor2[0]=curColor2[1]=curColor2[2]=0;
			for (i=y-halfWidth;i<y+halfWidth;i++){
				if ((i>=0)&&(i<height)){
					/*	Bottom */
					curColor[0]+=map[(x+i*width)*4+GlowR]*filter[fy];
					curColor[1]+=map[(x+i*width)*4+GlowG]*filter[fy];
					curColor[2]+=map[(x+i*width)*4+GlowB]*filter[fy];

					/*	Top */
					curColor2[0]+=map[(x+(height-1-i)*width) *
						4+GlowR]*filter[fy];
					curColor2[1]+=map[(x+(height-1-i)*width) *
						4+GlowG]*filter[fy];
					curColor2[2]+=map[(x+(height-1-i)*width) *
						4+GlowB]*filter[fy];
				}
				fy++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];
			temp[((x+(height-1-y)*width)*4)+GlowR]=curColor2[0];
			temp[((x+(height-1-y)*width)*4)+GlowG]=curColor2[1];
			temp[((x+(height-1-y)*width)*4)+GlowB]=curColor2[2];
		}
		/*	Do the main body */
		for (y=halfWidth;y<height-halfWidth;y++){
			index=(x+y*width)*4;
			fy=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			for (i=y-halfWidth;i<y+halfWidth;i++){
				curColor[0]+=map[(x+i*width)*4+GlowR]*filter[fy];
				curColor[1]+=map[(x+i*width)*4+GlowG]*filter[fy];
				curColor[2]+=map[(x+i*width)*4+GlowB]*filter[fy];
				fy++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];
		}
	}


	/*	Swap buffers */
	swap=temp;temp=map;map=swap;

	/*	Tidy up	*/
	MEM_freeN (filter);
	MEM_freeN (temp);
}

static void RVBlurBitmap2_float ( float* map, int width,int height,
				float blur, int quality)
/*	MUUUCCH better than the previous blur. */
/*	We do the blurring in two passes which is a whole lot faster. */
/*	I changed the math arount to implement an actual Gaussian */
/*	distribution. */
/*	*/
/*	Watch out though, it tends to misbehaven with large blur values on */
/*	a small bitmap.  Avoid avoid avoid. */
/*=============================== */
{
	float*	temp=NULL,*swap;
	float	*filter=NULL;
	int		x,y,i,fx,fy;
	int		index, ix, halfWidth;
	float	fval, k, curColor[3], curColor2[3], weight=0;

	/*	If we're not really blurring, bail out */
	if (blur<=0)
		return;

	/*	Allocate memory for the tempmap and the blur filter matrix */
	temp= MEM_mallocN( (width*height*4*sizeof(float)), "blurbitmaptemp");
	if (!temp)
		return;

	/*	Allocate memory for the filter elements */
	halfWidth = ((quality+1)*blur);
	filter = (float *)MEM_mallocN(sizeof(float)*halfWidth*2, "blurbitmapfilter");
	if (!filter){
		MEM_freeN (temp);
		return;
	}

	/*	Apparently we're calculating a bell curve */
	/*	based on the standard deviation (or radius) */
	/*	This code is based on an example */
	/*	posted to comp.graphics.algorithms by */
	/*	Blancmange (bmange@airdmhor.gen.nz) */

	k = -1.0f/(2.0f*(float)M_PI*blur*blur);

	for (ix = 0;ix< halfWidth;ix++){
		weight = (float)exp(k*(ix*ix));
		filter[halfWidth - ix] = weight;
		filter[halfWidth + ix] = weight;
	}
	filter[0] = weight;

	/*	Normalize the array */
	fval=0;
	for (ix = 0;ix< halfWidth*2;ix++)
		fval+=filter[ix];

	for (ix = 0;ix< halfWidth*2;ix++)
		filter[ix]/=fval;

	/*	Blur the rows */
	for (y=0;y<height;y++){
		/*	Do the left & right strips */
		for (x=0;x<halfWidth;x++){
			index=(x+y*width)*4;
			fx=0;
			curColor[0]=curColor[1]=curColor[2]=0.0f;
			curColor2[0]=curColor2[1]=curColor2[2]=0.0f;

			for (i=x-halfWidth;i<x+halfWidth;i++){
				if ((i>=0)&&(i<width)){
					curColor[0]+=map[(i+y*width)*4+GlowR]*filter[fx];
					curColor[1]+=map[(i+y*width)*4+GlowG]*filter[fx];
					curColor[2]+=map[(i+y*width)*4+GlowB]*filter[fx];

					curColor2[0]+=map[(width-1-i+y*width)*4+GlowR] *
						filter[fx];
					curColor2[1]+=map[(width-1-i+y*width)*4+GlowG] *
						filter[fx];
					curColor2[2]+=map[(width-1-i+y*width)*4+GlowB] *
						filter[fx];
				}
				fx++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];

			temp[((width-1-x+y*width)*4)+GlowR]=curColor2[0];
			temp[((width-1-x+y*width)*4)+GlowG]=curColor2[1];
			temp[((width-1-x+y*width)*4)+GlowB]=curColor2[2];

		}
		/*	Do the main body */
		for (x=halfWidth;x<width-halfWidth;x++){
			index=(x+y*width)*4;
			fx=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			for (i=x-halfWidth;i<x+halfWidth;i++){
				curColor[0]+=map[(i+y*width)*4+GlowR]*filter[fx];
				curColor[1]+=map[(i+y*width)*4+GlowG]*filter[fx];
				curColor[2]+=map[(i+y*width)*4+GlowB]*filter[fx];
				fx++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];
		}
	}

	/*	Swap buffers */
	swap=temp;temp=map;map=swap;


	/*	Blur the columns */
	for (x=0;x<width;x++){
		/*	Do the top & bottom strips */
		for (y=0;y<halfWidth;y++){
			index=(x+y*width)*4;
			fy=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			curColor2[0]=curColor2[1]=curColor2[2]=0;
			for (i=y-halfWidth;i<y+halfWidth;i++){
				if ((i>=0)&&(i<height)){
					/*	Bottom */
					curColor[0]+=map[(x+i*width)*4+GlowR]*filter[fy];
					curColor[1]+=map[(x+i*width)*4+GlowG]*filter[fy];
					curColor[2]+=map[(x+i*width)*4+GlowB]*filter[fy];

					/*	Top */
					curColor2[0]+=map[(x+(height-1-i)*width) *
						4+GlowR]*filter[fy];
					curColor2[1]+=map[(x+(height-1-i)*width) *
						4+GlowG]*filter[fy];
					curColor2[2]+=map[(x+(height-1-i)*width) *
						4+GlowB]*filter[fy];
				}
				fy++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];
			temp[((x+(height-1-y)*width)*4)+GlowR]=curColor2[0];
			temp[((x+(height-1-y)*width)*4)+GlowG]=curColor2[1];
			temp[((x+(height-1-y)*width)*4)+GlowB]=curColor2[2];
		}
		/*	Do the main body */
		for (y=halfWidth;y<height-halfWidth;y++){
			index=(x+y*width)*4;
			fy=0;
			curColor[0]=curColor[1]=curColor[2]=0;
			for (i=y-halfWidth;i<y+halfWidth;i++){
				curColor[0]+=map[(x+i*width)*4+GlowR]*filter[fy];
				curColor[1]+=map[(x+i*width)*4+GlowG]*filter[fy];
				curColor[2]+=map[(x+i*width)*4+GlowB]*filter[fy];
				fy++;
			}
			temp[index+GlowR]=curColor[0];
			temp[index+GlowG]=curColor[1];
			temp[index+GlowB]=curColor[2];
		}
	}


	/*	Swap buffers */
	swap=temp;temp=map;map=swap;

	/*	Tidy up	*/
	MEM_freeN (filter);
	MEM_freeN (temp);
}


/*	Adds two bitmaps and puts the results into a third map. */
/*	C must have been previously allocated but it may be A or B. */
/*	We clamp values to 255 to prevent weirdness */
/*=============================== */
static void RVAddBitmaps_byte (unsigned char* a, unsigned char* b,
				unsigned char* c, int width, int height)
{
	int	x,y,index;

	for (y=0;y<height;y++){
		for (x=0;x<width;x++){
			index=(x+y*width)*4;
			c[index+GlowR]=MIN2(255,a[index+GlowR]+b[index+GlowR]);
			c[index+GlowG]=MIN2(255,a[index+GlowG]+b[index+GlowG]);
			c[index+GlowB]=MIN2(255,a[index+GlowB]+b[index+GlowB]);
			c[index+GlowA]=MIN2(255,a[index+GlowA]+b[index+GlowA]);
		}
	}
}

static void RVAddBitmaps_float (float* a, float* b, float* c,
				int width, int height)
{
	int	x,y,index;

	for (y=0;y<height;y++){
		for (x=0;x<width;x++){
			index=(x+y*width)*4;
			c[index+GlowR]= MIN2(1.0f, a[index+GlowR]+b[index+GlowR]);
			c[index+GlowG]= MIN2(1.0f, a[index+GlowG]+b[index+GlowG]);
			c[index+GlowB]= MIN2(1.0f, a[index+GlowB]+b[index+GlowB]);
			c[index+GlowA]= MIN2(1.0f, a[index+GlowA]+b[index+GlowA]);
		}
	}
}

/*	For each pixel whose total luminance exceeds the threshold, */
/*	Multiply it's value by BOOST and add it to the output map */
static void RVIsolateHighlights_byte (unsigned char* in, unsigned char* out,
					int width, int height, int threshold,
					float boost, float clamp)
{
	int	x,y,index;
	int	intensity;


	for(y=0;y< height;y++) {
		for (x=0;x< width;x++) {
			index= (x+y*width)*4;

			/*	Isolate the intensity */
			intensity=(in[index+GlowR]+in[index+GlowG]+in[index+GlowB]-threshold);
			if (intensity>0){
				out[index+GlowR]=MIN2(255*clamp, (in[index+GlowR]*boost*intensity)/255);
				out[index+GlowG]=MIN2(255*clamp, (in[index+GlowG]*boost*intensity)/255);
				out[index+GlowB]=MIN2(255*clamp, (in[index+GlowB]*boost*intensity)/255);
				out[index+GlowA]=MIN2(255*clamp, (in[index+GlowA]*boost*intensity)/255);
			} else{
				out[index+GlowR]=0;
				out[index+GlowG]=0;
				out[index+GlowB]=0;
				out[index+GlowA]=0;
			}
		}
	}
}

static void RVIsolateHighlights_float (float* in, float* out,
					int width, int height, float threshold,
					float boost, float clamp)
{
	int		x,y,index;
	float	intensity;


	for(y=0;y< height;y++) {
		for (x=0;x< width;x++) {
			index= (x+y*width)*4;

			/*	Isolate the intensity */
			intensity=(in[index+GlowR]+in[index+GlowG]+in[index+GlowB]-threshold);
			if (intensity>0){
				out[index+GlowR]=MIN2(clamp, (in[index+GlowR]*boost*intensity));
				out[index+GlowG]=MIN2(clamp, (in[index+GlowG]*boost*intensity));
				out[index+GlowB]=MIN2(clamp, (in[index+GlowB]*boost*intensity));
				out[index+GlowA]=MIN2(clamp, (in[index+GlowA]*boost*intensity));
			} else{
				out[index+GlowR]=0;
				out[index+GlowG]=0;
				out[index+GlowB]=0;
				out[index+GlowA]=0;
			}
		}
	}
}

static void init_glow_effect(Sequence *seq)
{
	GlowVars *glow;

	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = MEM_callocN(sizeof(struct GlowVars), "glowvars");

	glow = (GlowVars *)seq->effectdata;
	glow->fMini = 0.25;
	glow->fClamp = 1.0;
	glow->fBoost = 0.5;
	glow->dDist = 3.0;
	glow->dQuality = 3;
	glow->bNoComp = 0;
}

static int num_inputs_glow(void)
{
	return 1;
}

static void free_glow_effect(Sequence *seq)
{
	if(seq->effectdata)MEM_freeN(seq->effectdata);
	seq->effectdata = NULL;
}

static void copy_glow_effect(Sequence *dst, Sequence *src)
{
	dst->effectdata = MEM_dupallocN(src->effectdata);
}

static void do_glow_effect_byte(Sequence *seq, int render_size, float facf0, float UNUSED(facf1),
				int x, int y, char *rect1, char *UNUSED(rect2), char *out)
{
	unsigned char *outbuf=(unsigned char *)out;
	unsigned char *inbuf=(unsigned char *)rect1;
	GlowVars *glow = (GlowVars *)seq->effectdata;

	RVIsolateHighlights_byte(inbuf, outbuf , x, y, glow->fMini*765, glow->fBoost * facf0, glow->fClamp);
	RVBlurBitmap2_byte (outbuf, x, y, glow->dDist * (render_size / 100.0f),glow->dQuality);
	if (!glow->bNoComp)
		RVAddBitmaps_byte (inbuf , outbuf, outbuf, x, y);
}

static void do_glow_effect_float(Sequence *seq, int render_size, float facf0, float UNUSED(facf1),
				int x, int y, float *rect1, float *UNUSED(rect2), float *out)
{
	float *outbuf = out;
	float *inbuf = rect1;
	GlowVars *glow = (GlowVars *)seq->effectdata;

	RVIsolateHighlights_float(inbuf, outbuf , x, y, glow->fMini*3.0f, glow->fBoost * facf0, glow->fClamp);
	RVBlurBitmap2_float (outbuf, x, y, glow->dDist * (render_size / 100.0f),glow->dQuality);
	if (!glow->bNoComp)
		RVAddBitmaps_float (inbuf , outbuf, outbuf, x, y);
}

static struct ImBuf * do_glow_effect(
	SeqRenderData context, Sequence *seq, float UNUSED(cfra),
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);

	int render_size = 100*context.rectx/context.scene->r.xsch;

	if (out->rect_float) {
		do_glow_effect_float(seq, render_size,
				facf0, facf1,
				context.rectx, context.recty,
				ibuf1->rect_float, ibuf2->rect_float,
				out->rect_float);
	} else {
		do_glow_effect_byte(seq, render_size,
				facf0, facf1,
				context.rectx, context.recty,
				(char*) ibuf1->rect, (char*) ibuf2->rect,
				(char*) out->rect);
	}

	return out;
}

/* setup */
void SeqConfigHandle_glow(struct SeqEffectHandle *hndl)
{
	hndl->init = init_glow_effect;
	hndl->num_inputs = num_inputs_glow;
	hndl->free = free_glow_effect;
	hndl->copy = copy_glow_effect;
	hndl->execute = do_glow_effect;
}


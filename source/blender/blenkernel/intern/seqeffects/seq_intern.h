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

/** \file blender/blenkernel/intern/seqeffects/seq_intern.h
 *  \ingroup bke
 */


/* from seqeffects.c */
extern struct ImBuf * prepare_effect_imbufs(SeqRenderData context, struct ImBuf *ibuf1, struct ImBuf *ibuf2, struct ImBuf *ibuf3);


/* The following are used by more than one effect - may or may not be useful for others */
/* from seq_alpha.c */
extern void do_alphaover_effect_byte(float facf0, float facf1, int x, int y, char *rect1, char *rect2, char *out);
extern void do_alphaover_effect_float(float facf0, float facf1, int x, int y, float *rect1, float *rect2, float *out);
extern void do_alphaunder_effect_byte(float facf0, float facf1, int x, int y, char *rect1, char *rect2, char *out);
extern void do_alphaunder_effect_float(float facf0, float facf1, int x, int y, float *rect1, float *rect2, float *out);

/* from seq_cross.c */
extern struct ImBuf* do_cross_effect(SeqRenderData context, Sequence *UNUSED(seq),
			float UNUSED(cfra), float facf0, float facf1, struct ImBuf *ibuf1,
			struct ImBuf *ibuf2, struct ImBuf *ibuf3);


/*
 * $Id$
 *
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. The Blender
 * Foundation also sells licenses for use in proprietary software under
 * the Blender License.  See http://www.blender.org/BL/ for information
 * about this.
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
 * The Original Code is Copyright (C) 2004 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Blender Foundation (2008).
 * Shane Ambler 2011
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef SEQ_EFFECTS_H
#define SEQ_EFFECTS_H

/** \file SEQ_effects.h
 *  \ingroup seq
 */

struct ImBuf;
struct Scene;
struct Sequence;
struct SeqRenderData;

/* Wipe effect */
enum {DO_SINGLE_WIPE, DO_DOUBLE_WIPE, DO_BOX_WIPE, DO_CROSS_WIPE,
	DO_IRIS_WIPE,DO_CLOCK_WIPE};


struct SeqEffectHandle {
	/* constructors & destructor */
	/* init & init_plugin are _only_ called on first creation */
	void (*init)(struct Sequence *seq);
	void (*init_plugin)(struct Sequence *seq, const char *fname);

	/* number of input strips needed
		(called directly after construction) */
	int (*num_inputs)(void);

	/* load is called first time after readblenfile in
		get_sequence_effect automatically */
	void (*load)(struct Sequence *seq);

	/* duplicate */
	void (*copy)(struct Sequence *dst, struct Sequence *src);

	/* destruct */
	void (*free)(struct Sequence *seq);

	/* returns: -1: no input needed,
	0: no early out,
	1: out = ibuf1,
	2: out = ibuf2 */
	int (*early_out)(struct Sequence *seq, float facf0, float facf1);

	/* stores the y-range of the effect IPO */
	void (*store_icu_yrange)(struct Sequence * seq,
                                 short adrcode, float *ymin, float *ymax);

	/* stores the default facf0 and facf1 if no IPO is present */
	void (*get_default_fac)(struct Sequence *seq, float cfra,
                                float * facf0, float * facf1);

	/* execute the effect
           sequence effects are only required to either support
           float-rects or byte-rects
           (mixed cases are handled one layer up...) */

	struct ImBuf* (*execute)(
		SeqRenderData context,
		struct Sequence *seq, float cfra,
		float facf0, float facf1,
		struct ImBuf *ibuf1, struct ImBuf *ibuf2,
		struct ImBuf *ibuf3);
};

/* **********************************************************************
   seqeffects.c

   Sequencer effect strip managment functions
   **********************************************************************
*/

struct SeqEffectHandle get_sequence_effect(struct Sequence *seq);
struct SeqEffectHandle get_sequence_blend(struct Sequence *seq);
int get_sequence_effect_num_inputs(int seq_type);
void sequence_effect_speed_rebuild_map(struct Scene *scene, struct Sequence *seq, int force);

#endif // SEQ_EFFECTS_H

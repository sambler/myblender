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

/** \file blender/blenkernel/intern/seqeffects/seq_multicam.c
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
   MULTICAM
   ********************************************************************** */

/* no effect inputs for multicam, we use give_ibuf_seq */
static int num_inputs_multicam(void)
{
	return 0;
}

static int early_out_multicam(struct Sequence *UNUSED(seq), float UNUSED(facf0), float UNUSED(facf1))
{
	return -1;
}

static struct ImBuf * do_multicam(
	SeqRenderData context, Sequence *seq, float cfra,
	float UNUSED(facf0), float UNUSED(facf1),
	struct ImBuf *UNUSED(ibuf1), struct ImBuf *UNUSED(ibuf2),
	struct ImBuf *UNUSED(ibuf3))
{
	struct ImBuf * i;
	struct ImBuf * out;
	Editing * ed;
	ListBase * seqbasep;

	if (seq->multicam_source == 0 || seq->multicam_source >= seq->machine) {
		return NULL;
	}

	ed = context.scene->ed;
	if (!ed) {
		return NULL;
	}
	seqbasep = seq_seqbase(&ed->seqbase, seq);
	if (!seqbasep) {
		return NULL;
	}

	i = give_ibuf_seqbase(context, cfra, seq->multicam_source, seqbasep);
	if (!i) {
		return NULL;
	}

	if (input_have_to_preprocess(context, seq, cfra)) {
		out = IMB_dupImBuf(i);
		IMB_freeImBuf(i);
	} else {
		out = i;
	}

	return out;
}

/* setup */
void SeqConfigHandle_multicam(struct SeqEffectHandle *hndl)
{
	hndl->num_inputs = num_inputs_multicam;
	hndl->early_out = early_out_multicam;
	hndl->execute = do_multicam;
}


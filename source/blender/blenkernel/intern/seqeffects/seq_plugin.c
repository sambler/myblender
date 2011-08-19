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

/** \file blender/blenkernel/intern/seqeffects/seq_plugin.c
 *  \ingroup bke
 */

#include <stdlib.h>

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"

#include "BKE_plugin_types.h"
#include "BKE_sequencer.h"
#include "BKE_utildefines.h"

#include "IMB_imbuf_types.h"
#include "IMB_imbuf.h"

#include "seq_intern.h"

#define INT	96
#define FLO	128

/* TODO: handle sequencer plugin errors */
static void error(const char *UNUSED(error), ...) {}

/* **********************************************************************
   PLUGINS
   ********************************************************************** */

static void open_plugin_seq(PluginSeq *pis, const char *seqname)
{
	int (*version)(void);
	void* (*alloc_private)(void);
	char *cp;

	/* to be sure: (is tested for) */
	pis->doit= NULL;
	pis->pname= NULL;
	pis->varstr= NULL;
	pis->cfra= NULL;
	pis->version= 0;
	pis->instance_private_data = NULL;

	/* clear the error list */
	BLI_dynlib_get_error_as_string(NULL);

	/* if(pis->handle) BLI_dynlib_close(pis->handle); */
	/* pis->handle= 0; */

	/* open the needed object */
	pis->handle= BLI_dynlib_open(pis->name);
	if(test_dlerr(pis->name, pis->name)) return;

	if (pis->handle != NULL) {
		/* find the address of the version function */
		version= (int (*)(void))BLI_dynlib_find_symbol(pis->handle, "plugin_seq_getversion");
		if (test_dlerr(pis->name, "plugin_seq_getversion")) return;

		if (version != NULL) {
			pis->version= version();
			if (pis->version >= 2 && pis->version <= 6) {
				int (*info_func)(PluginInfo *);
				PluginInfo *info= (PluginInfo*) MEM_mallocN(sizeof(PluginInfo), "plugin_info");

				info_func= (int (*)(PluginInfo *))BLI_dynlib_find_symbol(pis->handle, "plugin_getinfo");

				if(info_func == NULL) error("No info func");
				else {
					info_func(info);

					pis->pname= info->name;
					pis->vars= info->nvars;
					pis->cfra= info->cfra;

					pis->varstr= info->varstr;

					pis->doit= (void(*)(void))info->seq_doit;
					if (info->init)
						info->init();
				}
				MEM_freeN(info);

				cp= BLI_dynlib_find_symbol(pis->handle, "seqname");
				if(cp) strncpy(cp, seqname, 21);
			} else {
				printf ("Plugin returned unrecognized version number\n");
				return;
			}
		}
		alloc_private = (void* (*)(void))BLI_dynlib_find_symbol(
			pis->handle, "plugin_seq_alloc_private_data");
		if (alloc_private) {
			pis->instance_private_data = alloc_private();
		}

		pis->current_private_data = (void**)
			BLI_dynlib_find_symbol(
				pis->handle, "plugin_private_data");
	}
}

static PluginSeq *add_plugin_seq(const char *str, const char *seqname)
{
	PluginSeq *pis;
	VarStruct *varstr;
	int a;

	pis= MEM_callocN(sizeof(PluginSeq), "PluginSeq");

	strncpy(pis->name, str, FILE_MAXDIR+FILE_MAXFILE);
	open_plugin_seq(pis, seqname);

	if(pis->doit==NULL) {
		if(pis->handle==NULL) error("no plugin: %s", str);
		else error("in plugin: %s", str);
		MEM_freeN(pis);
		return NULL;
	}

	/* default values */
	varstr= pis->varstr;
	for(a=0; a<pis->vars; a++, varstr++) {
		if( (varstr->type & FLO)==FLO)
			pis->data[a]= varstr->def;
		else if( (varstr->type & INT)==INT)
			*((int *)(pis->data+a))= (int) varstr->def;
	}

	return pis;
}

static void free_plugin_seq(PluginSeq *pis)
{
	if(pis==NULL) return;

	/* no BLI_dynlib_close: same plugin can be opened multiple times with 1 handle */

	if (pis->instance_private_data) {
		void (*free_private)(void *);

		free_private = (void (*)(void *))BLI_dynlib_find_symbol(
			pis->handle, "plugin_seq_free_private_data");
		if (free_private) {
			free_private(pis->instance_private_data);
		}
	}

	MEM_freeN(pis);
}

static void init_plugin(Sequence * seq, const char * fname)
{
	seq->plugin= (PluginSeq *)add_plugin_seq(fname, seq->name+2);
}

/*
 * FIXME: should query plugin! Could be generator, that needs zero inputs...
 */
static int num_inputs_plugin(void)
{
	return 1;
}

static void load_plugin(Sequence * seq)
{
	if (seq) {
		open_plugin_seq(seq->plugin, seq->name+2);
	}
}

static void copy_plugin(Sequence * dst, Sequence * src)
{
	if(src->plugin) {
		dst->plugin= MEM_dupallocN(src->plugin);
		open_plugin_seq(dst->plugin, dst->name+2);
	}
}

static ImBuf * IMB_cast_away_list(ImBuf * i)
{
	if (!i) {
		return NULL;
	}
	return (ImBuf*) (((void**) i) + 2);
}

static struct ImBuf * do_plugin_effect(
	SeqRenderData context, Sequence *seq, float cfra,
	float facf0, float facf1,
	struct ImBuf *ibuf1, struct ImBuf *ibuf2,
	struct ImBuf *ibuf3)
{
	char *cp;
	int float_rendering;
	int use_temp_bufs = 0; /* Are needed since blur.c (and maybe some other
				old plugins) do very bad stuff
				with imbuf-internals */

	struct ImBuf * out = prepare_effect_imbufs(context,ibuf1, ibuf2, ibuf3);
	int x = context.rectx;
	int y = context.recty;

	if(seq->plugin && seq->plugin->doit) {

		if(seq->plugin->cfra)
			*(seq->plugin->cfra)= cfra;

		cp = BLI_dynlib_find_symbol(
			seq->plugin->handle, "seqname");

		if(cp) strncpy(cp, seq->name+2, 22);

		if (seq->plugin->current_private_data) {
			*seq->plugin->current_private_data
				= seq->plugin->instance_private_data;
		}

		float_rendering = (out->rect_float != NULL);

		if (seq->plugin->version<=3 && float_rendering) {
			use_temp_bufs = 1;

			if (ibuf1) {
				ibuf1 = IMB_dupImBuf(ibuf1);
				IMB_rect_from_float(ibuf1);
				imb_freerectfloatImBuf(ibuf1);
				ibuf1->flags &= ~IB_rectfloat;
			}
			if (ibuf2) {
				ibuf2 = IMB_dupImBuf(ibuf2);
				IMB_rect_from_float(ibuf2);
				imb_freerectfloatImBuf(ibuf2);
				ibuf2->flags &= ~IB_rectfloat;
			}
			if (ibuf3) {
				ibuf3 = IMB_dupImBuf(ibuf3);
				IMB_rect_from_float(ibuf3);
				imb_freerectfloatImBuf(ibuf3);
				ibuf3->flags &= ~IB_rectfloat;
			}
			if (!out->rect) imb_addrectImBuf(out);
			imb_freerectfloatImBuf(out);
			out->flags &= ~IB_rectfloat;
		}

		if (seq->plugin->version<=2) {
			if(ibuf1) IMB_convert_rgba_to_abgr(ibuf1);
			if(ibuf2) IMB_convert_rgba_to_abgr(ibuf2);
			if(ibuf3) IMB_convert_rgba_to_abgr(ibuf3);
		}

		if (seq->plugin->version<=4) {
			((SeqDoit)seq->plugin->doit)(
				seq->plugin->data, facf0, facf1, x, y,
				IMB_cast_away_list(ibuf1),
				IMB_cast_away_list(ibuf2),
				IMB_cast_away_list(out),
				IMB_cast_away_list(ibuf3));
		} else {
			((SeqDoit)seq->plugin->doit)(
				seq->plugin->data, facf0, facf1, x, y,
				ibuf1, ibuf2, out, ibuf3);
		}

		if (seq->plugin->version<=2) {
			if (!use_temp_bufs) {
				if(ibuf1) IMB_convert_rgba_to_abgr(ibuf1);
				if(ibuf2) IMB_convert_rgba_to_abgr(ibuf2);
				if(ibuf3) IMB_convert_rgba_to_abgr(ibuf3);
			}
			IMB_convert_rgba_to_abgr(out);
		}
		if (seq->plugin->version<=3 && float_rendering) {
			IMB_float_from_rect_simple(out);
		}

		if (use_temp_bufs) {
			if (ibuf1) IMB_freeImBuf(ibuf1);
			if (ibuf2) IMB_freeImBuf(ibuf2);
			if (ibuf3) IMB_freeImBuf(ibuf3);
		}
	}
	return out;
}

static int do_plugin_early_out(struct Sequence *UNUSED(seq),
				float UNUSED(facf0), float UNUSED(facf1))
{
	return 0;
}

static void free_plugin(struct Sequence * seq)
{
	free_plugin_seq(seq->plugin);
	seq->plugin = NULL;
}

/* setup */
void SeqConfigHandle_plugin(struct SeqEffectHandle *hndl)
{
	hndl->init_plugin = init_plugin;
	hndl->num_inputs = num_inputs_plugin;
	hndl->load = load_plugin;
	hndl->free = free_plugin;
	hndl->copy = copy_plugin;
	hndl->execute = do_plugin_effect;
	hndl->early_out = do_plugin_early_out;
}


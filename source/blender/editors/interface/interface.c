/*
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
 * Contributor(s): Blender Foundation 2002-2008, full recode.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/interface/interface.c
 *  \ingroup edinterface
 */


#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stddef.h>  /* offsetof() */

#include "MEM_guardedalloc.h"

#include "DNA_scene_types.h"
#include "DNA_screen_types.h"
#include "DNA_userdef_types.h"

#include "BLI_math.h"
#include "BLI_listbase.h"
#include "BLI_string.h"
#include "BLI_string_utf8.h"
#include "BLI_path_util.h"
#include "BLI_rect.h"

#include "BLI_dynstr.h"
#include "BLI_utildefines.h"

#include "BKE_context.h"
#include "BKE_library.h"
#include "BKE_unit.h"
#include "BKE_screen.h"
#include "BKE_idprop.h"

#include "BIF_gl.h"

#include "BLF_api.h"
#include "BLF_translation.h"

#include "UI_interface.h"

#include "IMB_imbuf.h"

#include "WM_api.h"
#include "WM_types.h"
#include "wm_subwindow.h"
#include "wm_window.h"

#include "RNA_access.h"

#include "BPY_extern.h"

#include "IMB_colormanagement.h"

#include "interface_intern.h"

#define PRECISION_FLOAT_MAX 6
#define PRECISION_FLOAT_MAX_POW 1000000 /* pow(10, PRECISION_FLOAT_MAX)  */

/* avoid unneeded calls to ui_get_but_val */
#define UI_BUT_VALUE_UNSET DBL_MAX
#define UI_GET_BUT_VALUE_INIT(_but, _value) if (_value == DBL_MAX) {  (_value) = ui_get_but_val(_but); } (void)0

/* 
 * a full doc with API notes can be found in bf-blender/trunk/blender/doc/guides/interface_API.txt
 * 
 * uiBlahBlah()		external function
 * ui_blah_blah()	internal function
 */

static void ui_free_but(const bContext *C, uiBut *but);

/* ************* window matrix ************** */

void ui_block_to_window_fl(const ARegion *ar, uiBlock *block, float *x, float *y)
{
	float gx, gy;
	int sx, sy, getsizex, getsizey;

	getsizex = BLI_rcti_size_x(&ar->winrct) + 1;
	getsizey = BLI_rcti_size_y(&ar->winrct) + 1;
	sx = ar->winrct.xmin;
	sy = ar->winrct.ymin;

	gx = *x;
	gy = *y;

	if (block->panel) {
		gx += block->panel->ofsx;
		gy += block->panel->ofsy;
	}

	*x = ((float)sx) + ((float)getsizex) * (0.5f + 0.5f * (gx * block->winmat[0][0] + gy * block->winmat[1][0] + block->winmat[3][0]));
	*y = ((float)sy) + ((float)getsizey) * (0.5f + 0.5f * (gx * block->winmat[0][1] + gy * block->winmat[1][1] + block->winmat[3][1]));
}

void ui_block_to_window(const ARegion *ar, uiBlock *block, int *x, int *y)
{
	float fx, fy;

	fx = *x;
	fy = *y;

	ui_block_to_window_fl(ar, block, &fx, &fy);

	*x = (int)(fx + 0.5f);
	*y = (int)(fy + 0.5f);
}

void ui_block_to_window_rct(const ARegion *ar, uiBlock *block, const rctf *graph, rcti *winr)
{
	rctf tmpr;

	tmpr = *graph;
	ui_block_to_window_fl(ar, block, &tmpr.xmin, &tmpr.ymin);
	ui_block_to_window_fl(ar, block, &tmpr.xmax, &tmpr.ymax);

	BLI_rcti_rctf_copy(winr, &tmpr);
}

void ui_window_to_block_fl(const ARegion *ar, uiBlock *block, float *x, float *y)   /* for mouse cursor */
{
	float a, b, c, d, e, f, px, py;
	int sx, sy, getsizex, getsizey;

	getsizex = BLI_rcti_size_x(&ar->winrct) + 1;
	getsizey = BLI_rcti_size_y(&ar->winrct) + 1;
	sx = ar->winrct.xmin;
	sy = ar->winrct.ymin;

	a = 0.5f * ((float)getsizex) * block->winmat[0][0];
	b = 0.5f * ((float)getsizex) * block->winmat[1][0];
	c = 0.5f * ((float)getsizex) * (1.0f + block->winmat[3][0]);

	d = 0.5f * ((float)getsizey) * block->winmat[0][1];
	e = 0.5f * ((float)getsizey) * block->winmat[1][1];
	f = 0.5f * ((float)getsizey) * (1.0f + block->winmat[3][1]);

	px = *x - sx;
	py = *y - sy;

	*y =  (a * (py - f) + d * (c - px)) / (a * e - d * b);
	*x = (px - b * (*y) - c) / a;

	if (block->panel) {
		*x -= block->panel->ofsx;
		*y -= block->panel->ofsy;
	}
}

void ui_window_to_block(const ARegion *ar, uiBlock *block, int *x, int *y)
{
	float fx, fy;

	fx = *x;
	fy = *y;

	ui_window_to_block_fl(ar, block, &fx, &fy);

	*x = (int)(fx + 0.5f);
	*y = (int)(fy + 0.5f);
}

void ui_window_to_region(const ARegion *ar, int *x, int *y)
{
	*x -= ar->winrct.xmin;
	*y -= ar->winrct.ymin;
}

/* ******************* block calc ************************* */

void ui_block_translate(uiBlock *block, int x, int y)
{
	uiBut *but;

	for (but = block->buttons.first; but; but = but->next) {
		BLI_rctf_translate(&but->rect, x, y);
	}

	BLI_rctf_translate(&block->rect, x, y);
}

static void ui_text_bounds_block(uiBlock *block, float offset)
{
	uiStyle *style = UI_GetStyle();
	uiBut *bt;
	int i = 0, j, x1addval = offset, nextcol;
	int lastcol = 0, col = 0;
	
	uiStyleFontSet(&style->widget);
	
	for (bt = block->buttons.first; bt; bt = bt->next) {
		if (bt->type != SEPR) {
			j = BLF_width(style->widget.uifont_id, bt->drawstr);

			if (j > i) i = j;
		}

		if (bt->next && bt->rect.xmin < bt->next->rect.xmin)
			lastcol++;
	}

	/* cope with multi collumns */
	bt = block->buttons.first;
	while (bt) {
		if (bt->next && bt->rect.xmin < bt->next->rect.xmin) {
			nextcol = 1;
			col++;
		}
		else nextcol = 0;
		
		bt->rect.xmin = x1addval;
		bt->rect.xmax = bt->rect.xmin + i + block->bounds;
		
		if (col == lastcol) {
			bt->rect.xmax = maxf(bt->rect.xmax, offset + block->minbounds);
		}

		ui_check_but(bt);  /* clips text again */
		
		if (nextcol)
			x1addval += i + block->bounds;
		
		bt = bt->next;
	}
}

void ui_bounds_block(uiBlock *block)
{
	uiBut *bt;
	int xof;
	
	if (block->buttons.first == NULL) {
		if (block->panel) {
			block->rect.xmin = 0.0; block->rect.xmax = block->panel->sizex;
			block->rect.ymin = 0.0; block->rect.ymax = block->panel->sizey;
		}
	}
	else {
	
		BLI_rctf_init_minmax(&block->rect);

		for (bt = block->buttons.first; bt; bt = bt->next) {
			BLI_rctf_union(&block->rect, &bt->rect);
		}
		
		block->rect.xmin -= block->bounds;
		block->rect.ymin -= block->bounds;
		block->rect.xmax += block->bounds;
		block->rect.ymax += block->bounds;
	}

	block->rect.xmax = block->rect.xmin + maxf(BLI_rctf_size_x(&block->rect), block->minbounds);

	/* hardcoded exception... but that one is annoying with larger safety */ 
	bt = block->buttons.first;
	if (bt && strncmp(bt->str, "ERROR", 5) == 0) xof = 10;
	else xof = 40;

	block->safety.xmin = block->rect.xmin - xof;
	block->safety.ymin = block->rect.ymin - xof;
	block->safety.xmax = block->rect.xmax + xof;
	block->safety.ymax = block->rect.ymax + xof;
}

static void ui_centered_bounds_block(const bContext *C, uiBlock *block)
{
	wmWindow *window = CTX_wm_window(C);
	int xmax, ymax;
	int startx, starty;
	int width, height;
	
	/* note: this is used for the splash where window bounds event has not been
	 * updated by ghost, get the window bounds from ghost directly */

	// wm_window_get_size(window, &xmax, &ymax);
	wm_window_get_size_ghost(window, &xmax, &ymax);
	
	ui_bounds_block(block);
	
	width  = BLI_rctf_size_x(&block->rect);
	height = BLI_rctf_size_y(&block->rect);
	
	startx = (xmax * 0.5f) - (width * 0.5f);
	starty = (ymax * 0.5f) - (height * 0.5f);
	
	ui_block_translate(block, startx - block->rect.xmin, starty - block->rect.ymin);
	
	/* now recompute bounds and safety */
	ui_bounds_block(block);
	
}
static void ui_popup_bounds_block(const bContext *C, uiBlock *block, eBlockBoundsCalc bounds_calc)
{
	wmWindow *window = CTX_wm_window(C);
	int startx, starty, endx, endy, width, height, oldwidth, oldheight;
	int oldbounds, xmax, ymax;

	oldbounds = block->bounds;

	/* compute mouse position with user defined offset */
	ui_bounds_block(block);
	
	wm_window_get_size(window, &xmax, &ymax);

	oldwidth  = BLI_rctf_size_x(&block->rect);
	oldheight = BLI_rctf_size_y(&block->rect);

	/* first we ensure wide enough text bounds */
	if (bounds_calc == UI_BLOCK_BOUNDS_POPUP_MENU) {
		if (block->flag & UI_BLOCK_LOOP) {
			block->bounds = 50;
			ui_text_bounds_block(block, block->rect.xmin);
		}
	}

	/* next we recompute bounds */
	block->bounds = oldbounds;
	ui_bounds_block(block);

	/* and we adjust the position to fit within window */
	width  = BLI_rctf_size_x(&block->rect);
	height = BLI_rctf_size_y(&block->rect);

	/* avoid divide by zero below, caused by calling with no UI, but better not crash */
	oldwidth = oldwidth > 0 ? oldwidth : MAX2(1, width);
	oldheight = oldheight > 0 ? oldheight : MAX2(1, height);

	/* offset block based on mouse position, user offset is scaled
	 * along in case we resized the block in ui_text_bounds_block */
	startx = window->eventstate->x + block->rect.xmin + (block->mx * width) / oldwidth;
	starty = window->eventstate->y + block->rect.ymin + (block->my * height) / oldheight;

	if (startx < 10)
		startx = 10;
	if (starty < 10)
		starty = 10;

	endx = startx + width;
	endy = starty + height;

	if (endx > xmax) {
		endx = xmax - 10;
		startx = endx - width;
	}
	if (endy > ymax - 20) {
		endy = ymax - 20;
		starty = endy - height;
	}

	ui_block_translate(block, startx - block->rect.xmin, starty - block->rect.ymin);

	/* now recompute bounds and safety */
	ui_bounds_block(block);
}

/* used for various cases */
void uiBoundsBlock(uiBlock *block, int addval)
{
	if (block == NULL)
		return;
	
	block->bounds = addval;
	block->bounds_type = UI_BLOCK_BOUNDS;
}

/* used for pulldowns */
void uiTextBoundsBlock(uiBlock *block, int addval)
{
	block->bounds = addval;
	block->bounds_type = UI_BLOCK_BOUNDS_TEXT;
}

/* used for block popups */
void uiPopupBoundsBlock(uiBlock *block, int addval, int mx, int my)
{
	block->bounds = addval;
	block->bounds_type = UI_BLOCK_BOUNDS_POPUP_MOUSE;
	block->mx = mx;
	block->my = my;
}

/* used for menu popups */
void uiMenuPopupBoundsBlock(uiBlock *block, int addval, int mx, int my)
{
	block->bounds = addval;
	block->bounds_type = UI_BLOCK_BOUNDS_POPUP_MENU;
	block->mx = mx;
	block->my = my;
}

/* used for centered popups, i.e. splash */
void uiCenteredBoundsBlock(uiBlock *block, int addval)
{
	block->bounds = addval;
	block->bounds_type = UI_BLOCK_BOUNDS_POPUP_CENTER;
}

void uiExplicitBoundsBlock(uiBlock *block, int minx, int miny, int maxx, int maxy)
{
	block->rect.xmin = minx;
	block->rect.ymin = miny;
	block->rect.xmax = maxx;
	block->rect.ymax = maxy;
	block->bounds_type = UI_BLOCK_BOUNDS_NONE;
}

/* ************** LINK LINE DRAWING  ************* */

/* link line drawing is not part of buttons or theme.. so we stick with it here */

static int ui_but_float_precision(uiBut *but, double value)
{
	int prec;

	/* first check if prec is 0 and fallback to a simple default */
	if ((prec = (int)but->a2) == 0) {
		prec = (but->hardmax < 10.001f) ? 3 : 2;
	}

	/* check on the number of decimal places need to display
	 * the number, this is so 0.00001 is not displayed as 0.00,
	 * _but_, this is only for small values si 10.0001 will not get
	 * the same treatment */
	if (value != 0.0 && (value = ABS(value)) < 0.1) {
		int value_i = (int)((value * PRECISION_FLOAT_MAX_POW) + 0.5);
		if (value_i != 0) {
			const int prec_span = 3; /* show: 0.01001, 5 would allow 0.0100001 for eg. */
			int test_prec;
			int prec_min = -1;
			int dec_flag = 0;
			int i = PRECISION_FLOAT_MAX;
			while (i && value_i) {
				if (value_i % 10) {
					dec_flag |= 1 << i;
					prec_min = i;
				}
				value_i /= 10;
				i--;
			}

			/* even though its a small value, if the second last digit is not 0, use it */
			test_prec = prec_min;

			dec_flag = (dec_flag >> (prec_min + 1)) & ((1 << prec_span) - 1);

			while (dec_flag) {
				test_prec++;
				dec_flag = dec_flag >> 1;
			}

			if (test_prec > prec) {
				prec = test_prec;
			}
		}
	}

	CLAMP(prec, 1, PRECISION_FLOAT_MAX);

	return prec;
}

static void ui_draw_linkline(uiLinkLine *line, int highlightActiveLines)
{
	rcti rect;

	if (line->from == NULL || line->to == NULL) return;
	
	rect.xmin = BLI_rctf_cent_x(&line->from->rect);
	rect.ymin = BLI_rctf_cent_y(&line->from->rect);
	rect.xmax = BLI_rctf_cent_x(&line->to->rect);
	rect.ymax = BLI_rctf_cent_y(&line->to->rect);
	
	if (line->flag & UI_SELECT)
		glColor3ub(100, 100, 100);
	else if (highlightActiveLines && ((line->from->flag & UI_ACTIVE) || (line->to->flag & UI_ACTIVE)))
		UI_ThemeColor(TH_TEXT_HI);
	else 
		glColor3ub(0, 0, 0);

	ui_draw_link_bezier(&rect);
}

static void ui_draw_links(uiBlock *block)
{
	uiBut *but;
	uiLinkLine *line;

	/* Draw the inactive lines (lines with neither button being hovered over).
	 * As we go, remember if we see any active or selected lines. */
	int foundselectline = FALSE;
	int foundactiveline = FALSE;
	for (but = block->buttons.first; but; but = but->next) {
		if (but->type == LINK && but->link) {
			for (line = but->link->lines.first; line; line = line->next) {
				if (!(line->from->flag & UI_ACTIVE) && !(line->to->flag & UI_ACTIVE))
					ui_draw_linkline(line, 0);
				else
					foundactiveline = TRUE;

				if ((line->from->flag & UI_SELECT) || (line->to->flag & UI_SELECT))
					foundselectline = TRUE;
			}
		}
	}	

	/* Draw any active lines (lines with either button being hovered over).
	 * Do this last so they appear on top of inactive lines. */
	if (foundactiveline) {
		for (but = block->buttons.first; but; but = but->next) {
			if (but->type == LINK && but->link) {
				for (line = but->link->lines.first; line; line = line->next) {
					if ((line->from->flag & UI_ACTIVE) || (line->to->flag & UI_ACTIVE))
						ui_draw_linkline(line, !foundselectline);
				}
			}
		}	
	}
}

/* ************** BLOCK ENDING FUNCTION ************* */

/* NOTE: if but->poin is allocated memory for every defbut, things fail... */
static int ui_but_equals_old(uiBut *but, uiBut *oldbut)
{
	/* various properties are being compared here, hopefully sufficient
	 * to catch all cases, but it is simple to add more checks later */
	if (but->retval != oldbut->retval) return 0;
	if (but->rnapoin.data != oldbut->rnapoin.data) return 0;
	if (but->rnaprop != oldbut->rnaprop)
		if (but->rnaindex != oldbut->rnaindex) return 0;
	if (but->func != oldbut->func) return 0;
	if (but->funcN != oldbut->funcN) return 0;
	if (oldbut->func_arg1 != oldbut && but->func_arg1 != oldbut->func_arg1) return 0;
	if (oldbut->func_arg2 != oldbut && but->func_arg2 != oldbut->func_arg2) return 0;
	if (!but->funcN && ((but->poin != oldbut->poin && (uiBut *)oldbut->poin != oldbut) || but->pointype != oldbut->pointype)) return 0;
	if (but->optype != oldbut->optype) return 0;

	return 1;
}

/* oldbut is being inserted in new block, so we use the lines from new button, and replace button pointers */
static void ui_but_update_linklines(uiBlock *block, uiBut *oldbut, uiBut *newbut)
{
	uiLinkLine *line;
	uiBut *but;
	
	/* if active button is LINK */
	if (newbut->type == LINK && newbut->link) {
		
		SWAP(uiLink *, oldbut->link, newbut->link);
		
		for (line = oldbut->link->lines.first; line; line = line->next) {
			if (line->to == newbut)
				line->to = oldbut;
			if (line->from == newbut)
				line->from = oldbut;
		}
	}		
	
	/* check all other button links */
	for (but = block->buttons.first; but; but = but->next) {
		if (but != newbut && but->type == LINK && but->link) {
			for (line = but->link->lines.first; line; line = line->next) {
				if (line->to == newbut)
					line->to = oldbut;
				if (line->from == newbut)
					line->from = oldbut;
			}
		}
	}
}

static int ui_but_update_from_old_block(const bContext *C, uiBlock *block, uiBut **butpp)
{
	uiBlock *oldblock;
	uiBut *oldbut, *but = *butpp;
	int found = 0;

	oldblock = block->oldblock;
	if (!oldblock)
		return found;

	for (oldbut = oldblock->buttons.first; oldbut; oldbut = oldbut->next) {
		if (ui_but_equals_old(oldbut, but)) {
			if (oldbut->active) {
#if 0
//				but->flag= oldbut->flag;
#else
				/* exception! redalert flag can't be update from old button. 
				 * perhaps it should only copy specific flags rather than all. */
//				but->flag= (oldbut->flag & ~UI_BUT_REDALERT) | (but->flag & UI_BUT_REDALERT);
#endif
//				but->active= oldbut->active;
//				but->pos= oldbut->pos;
//				but->ofs= oldbut->ofs;
//				but->editstr= oldbut->editstr;
//				but->editval= oldbut->editval;
//				but->editvec= oldbut->editvec;
//				but->editcoba= oldbut->editcoba;
//				but->editcumap= oldbut->editcumap;
//				but->selsta= oldbut->selsta;
//				but->selend= oldbut->selend;
//				but->softmin= oldbut->softmin;
//				but->softmax= oldbut->softmax;
//				but->linkto[0]= oldbut->linkto[0];
//				but->linkto[1]= oldbut->linkto[1];
				found = 1;
//				oldbut->active= NULL;
			
				/* move button over from oldblock to new block */
				BLI_remlink(&oldblock->buttons, oldbut);
				BLI_insertlink(&block->buttons, but, oldbut);
				oldbut->block = block;
				*butpp = oldbut;
				
				/* still stuff needs to be copied */
				oldbut->rect = but->rect;
				oldbut->context = but->context; /* set by Layout */
				
				/* typically the same pointers, but not on undo/redo */
				/* XXX some menu buttons store button itself in but->poin. Ugly */
				if (oldbut->poin != (char *)oldbut) {
					SWAP(char *, oldbut->poin, but->poin);
					SWAP(void *, oldbut->func_argN, but->func_argN);
				}
				
				/* copy hardmin for list rows to prevent 'sticking' highlight to mouse position
				 * when scrolling without moving mouse (see [#28432]) */
				if (ELEM(oldbut->type, ROW, LISTROW))
					oldbut->hardmax = but->hardmax;
				
				ui_but_update_linklines(block, oldbut, but);
				
				BLI_remlink(&block->buttons, but);
				ui_free_but(C, but);
				
				/* note: if layout hasn't been applied yet, it uses old button pointers... */
			}
			else {
				/* ensures one button can get activated, and in case the buttons
				 * draw are the same this gives O(1) lookup for each button */
				BLI_remlink(&oldblock->buttons, oldbut);
				ui_free_but(C, oldbut);
			}
			
			break;
		}
	}
	
	return found;
}

/* needed for temporarily rename buttons, such as in outliner or file-select,
 * they should keep calling uiDefButs to keep them alive */
/* returns 0 when button removed */
int uiButActiveOnly(const bContext *C, uiBlock *block, uiBut *but)
{
	uiBlock *oldblock;
	uiBut *oldbut;
	int activate = FALSE, found = FALSE, isactive = FALSE;
	
	oldblock = block->oldblock;
	if (!oldblock)
		activate = TRUE;
	else {
		for (oldbut = oldblock->buttons.first; oldbut; oldbut = oldbut->next) {
			if (ui_but_equals_old(oldbut, but)) {
				found = TRUE;
				
				if (oldbut->active)
					isactive = TRUE;
				
				break;
			}
		}
	}
	if ((activate == TRUE) || (found == FALSE)) {
		ui_button_activate_do((bContext *)C, CTX_wm_region(C), but);
	}
	else if ((found == TRUE) && (isactive == FALSE)) {
		BLI_remlink(&block->buttons, but);
		ui_free_but(C, but);
		return 0;
	}
	
	return 1;
}

/* use to check if we need to disable undo, but don't make any changes
 * returns FALSE if undo needs to be disabled. */
static int ui_but_is_rna_undo(uiBut *but)
{
	if (but->rnapoin.id.data) {
		/* avoid undo push for buttons who's ID are screen or wm level
		 * we could disable undo for buttons with no ID too but may have
		 * unforeseen consequences, so best check for ID's we _know_ are not
		 * handled by undo - campbell */
		ID *id = but->rnapoin.id.data;
		if (ID_CHECK_UNDO(id) == FALSE) {
			return FALSE;
		}
		else {
			return TRUE;
		}
	}
	else if (but->rnapoin.type && !RNA_struct_undo_check(but->rnapoin.type)) {
		return FALSE;
	}

	return TRUE;
}

/* assigns automatic keybindings to menu items for fast access
 * (underline key in menu) */
static void ui_menu_block_set_keyaccels(uiBlock *block)
{
	uiBut *but;

	unsigned int menu_key_mask = 0;
	unsigned char menu_key;
	const char *str_pt;
	int pass;
	int tot_missing = 0;

	/* only do it before bounding */
	if (block->rect.xmin != block->rect.xmax)
		return;

	for (pass = 0; pass < 2; pass++) {
		/* 2 Passes, on for first letter only, second for any letter if first fails
		 * fun first pass on all buttons so first word chars always get first priority */

		for (but = block->buttons.first; but; but = but->next) {
			if (!ELEM5(but->type, BUT, BUTM, MENU, BLOCK, PULLDOWN) || (but->flag & UI_HIDDEN)) {
				/* pass */
			}
			else if (but->menu_key == '\0') {
				if (but->str) {
					for (str_pt = but->str; *str_pt; ) {
						menu_key = tolower(*str_pt);
						if ((menu_key >= 'a' && menu_key <= 'z') && !(menu_key_mask & 1 << (menu_key - 'a'))) {
							menu_key_mask |= 1 << (menu_key - 'a');
							break;
						}

						if (pass == 0) {
							/* Skip to next delimiter on first pass (be picky) */
							while (isalpha(*str_pt))
								str_pt++;

							if (*str_pt)
								str_pt++;
						}
						else {
							/* just step over every char second pass and find first usable key */
							str_pt++;
						}
					}

					if (*str_pt) {
						but->menu_key = menu_key;
					}
					else {
						/* run second pass */
						tot_missing++;
					}

					/* if all keys have been used just exit, unlikely */
					if (menu_key_mask == (1 << 26) - 1) {
						return;
					}
				}
			}
		}

		/* check if second pass is needed */
		if (!tot_missing) {
			break;
		}
	}
}

/* XXX, this code will shorten any allocated string to 'UI_MAX_NAME_STR'
 * since this is really long its unlikely to be an issue,
 * but this could be supported */
void ui_but_add_shortcut(uiBut *but, const char *shortcut_str, const short do_strip)
{

	if (do_strip) {
		char *cpoin = strchr(but->str, '|');
		if (cpoin) {
			*cpoin = '\0';
		}
	}

	/* without this, just allow stripping of the shortcut */
	if (shortcut_str) {
		char *butstr_orig;

		if (but->str != but->strdata) {
			butstr_orig = but->str; /* free after using as source buffer */
		}
		else {
			butstr_orig = BLI_strdup(but->str);
		}
		BLI_snprintf(but->strdata,
		             sizeof(but->strdata),
		             "%s|%s",
		             butstr_orig, shortcut_str);
		MEM_freeN(butstr_orig);
		but->str = but->strdata;
		ui_check_but(but);
	}
}

static void ui_menu_block_set_keymaps(const bContext *C, uiBlock *block)
{
	uiBut *but;
	char buf[128];

	/* for menu's */
	MenuType *mt;
	IDProperty *prop_menu = NULL;
	IDProperty *prop_menu_name = NULL;

	/* only do it before bounding */
	if (block->rect.xmin != block->rect.xmax)
		return;

	for (but = block->buttons.first; but; but = but->next) {
		if (but->optype) {
			IDProperty *prop = (but->opptr) ? but->opptr->data : NULL;

			if (WM_key_event_operator_string(C, but->optype->idname, but->opcontext, prop, TRUE,
			                                 buf, sizeof(buf)))
			{
				ui_but_add_shortcut(but, buf, FALSE);
			}
		}
		else if ((mt = uiButGetMenuType(but))) {
			/* only allocate menu property once */
			if (prop_menu == NULL) {
				/* annoying, create a property */
				IDPropertyTemplate val = {0};
				prop_menu = IDP_New(IDP_GROUP, &val, __func__); /* dummy, name is unimportant  */
				IDP_AddToGroup(prop_menu, (prop_menu_name = IDP_NewString("", "name", sizeof(mt->idname))));
			}

			IDP_AssignString(prop_menu_name, mt->idname, sizeof(mt->idname));

			if (WM_key_event_operator_string(C, "WM_OT_call_menu", WM_OP_INVOKE_REGION_WIN, prop_menu, FALSE,
			                                 buf, sizeof(buf)))
			{
				ui_but_add_shortcut(but, buf, FALSE);
			}
		}
	}

	if (prop_menu) {
		IDP_FreeProperty(prop_menu);
		MEM_freeN(prop_menu);
	}

#undef UI_MENU_KEY_STR_CAT

}

void uiEndBlock(const bContext *C, uiBlock *block)
{
	uiBut *but;
	Scene *scene = CTX_data_scene(C);

	/* inherit flags from 'old' buttons that was drawn here previous, based
	 * on matching buttons, we need this to make button event handling non
	 * blocking, while still allowing buttons to be remade each redraw as it
	 * is expected by blender code */
	for (but = block->buttons.first; but; but = but->next) {
		if (ui_but_update_from_old_block(C, block, &but))
			ui_check_but(but);
		
		/* temp? Proper check for graying out */
		if (but->optype) {
			wmOperatorType *ot = but->optype;

			if (but->context)
				CTX_store_set((bContext *)C, but->context);

			if (ot == NULL || WM_operator_poll_context((bContext *)C, ot, but->opcontext) == 0) {
				but->flag |= UI_BUT_DISABLED;
				but->lock = TRUE;
			}

			if (but->context)
				CTX_store_set((bContext *)C, NULL);
		}

		ui_but_anim_flag(but, (scene) ? scene->r.cfra : 0.0f);
	}

	if (block->oldblock) {
		block->auto_open = block->oldblock->auto_open;
		block->auto_open_last = block->oldblock->auto_open_last;
		block->tooltipdisabled = block->oldblock->tooltipdisabled;

		block->oldblock = NULL;
	}

	/* handle pending stuff */
	if (block->layouts.first) {
		uiBlockLayoutResolve(block, NULL, NULL);
	}
	ui_block_do_align(block);
	if ((block->flag & UI_BLOCK_LOOP) && (block->flag & UI_BLOCK_NUMSELECT)) {
		ui_menu_block_set_keyaccels(block); /* could use a different flag to check */
	}

	if (block->flag & UI_BLOCK_LOOP) {
		ui_menu_block_set_keymaps(C, block);
	}
	
	/* after keymaps! */
	switch (block->bounds_type) {
		case UI_BLOCK_BOUNDS_NONE:
			break;
		case UI_BLOCK_BOUNDS:
			ui_bounds_block(block);
			break;
		case UI_BLOCK_BOUNDS_TEXT:
			ui_text_bounds_block(block, 0.0f);
			break;
		case UI_BLOCK_BOUNDS_POPUP_CENTER:
			ui_centered_bounds_block(C, block);
			break;

			/* fallback */
		case UI_BLOCK_BOUNDS_POPUP_MOUSE:
		case UI_BLOCK_BOUNDS_POPUP_MENU:
			ui_popup_bounds_block(C, block, block->bounds_type);
			break;
	}

	if (block->rect.xmin == 0.0f && block->rect.xmax == 0.0f) {
		uiBoundsBlock(block, 0);
	}
	if (block->flag & UI_BUT_ALIGN) {
		uiBlockEndAlign(block);
	}

	block->endblock = 1;
}

/* ************** BLOCK DRAWING FUNCTION ************* */

void ui_fontscale(short *points, float aspect)
{
	if (aspect < 0.9f || aspect > 1.1f) {
		float pointsf = *points;
		
		/* for some reason scaling fonts goes too fast compared to widget size */
		aspect = sqrt(aspect);
		pointsf /= aspect;
		
		if (aspect > 1.0f)
			*points = ceilf(pointsf);
		else
			*points = floorf(pointsf);
	}
}

/* project button or block (but==NULL) to pixels in regionspace */
static void ui_but_to_pixelrect(rcti *rect, const ARegion *ar, uiBlock *block, uiBut *but)
{
	float gx, gy;
	float getsizex, getsizey;
	
	getsizex = ar->winx;
	getsizey = ar->winy;

	gx = (but ? but->rect.xmin : block->rect.xmin) + (block->panel ? block->panel->ofsx : 0.0f);
	gy = (but ? but->rect.ymin : block->rect.ymin) + (block->panel ? block->panel->ofsy : 0.0f);
	
	rect->xmin = floorf(getsizex * (0.5f + 0.5f * (gx * block->winmat[0][0] + gy * block->winmat[1][0] + block->winmat[3][0])));
	rect->ymin = floorf(getsizey * (0.5f + 0.5f * (gx * block->winmat[0][1] + gy * block->winmat[1][1] + block->winmat[3][1])));
	
	gx = (but ? but->rect.xmax : block->rect.xmax) + (block->panel ? block->panel->ofsx : 0.0f);
	gy = (but ? but->rect.ymax : block->rect.ymax) + (block->panel ? block->panel->ofsy : 0.0f);
	
	rect->xmax = floorf(getsizex * (0.5f + 0.5f * (gx * block->winmat[0][0] + gy * block->winmat[1][0] + block->winmat[3][0])));
	rect->ymax = floorf(getsizey * (0.5f + 0.5f * (gx * block->winmat[0][1] + gy * block->winmat[1][1] + block->winmat[3][1])));

}

/* uses local copy of style, to scale things down, and allow widgets to change stuff */
void uiDrawBlock(const bContext *C, uiBlock *block)
{
	uiStyle style = *UI_GetStyle();  /* XXX pass on as arg */
	ARegion *ar;
	uiBut *but;
	rcti rect;
	int multisample_enabled;
	
	/* get menu region or area region */
	ar = CTX_wm_menu(C);
	if (!ar)
		ar = CTX_wm_region(C);

	if (!block->endblock)
		uiEndBlock(C, block);

	/* disable AA, makes widgets too blurry */
	multisample_enabled = glIsEnabled(GL_MULTISAMPLE_ARB);
	if (multisample_enabled)
		glDisable(GL_MULTISAMPLE_ARB);

	/* we set this only once */
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	/* scale fonts */
	ui_fontscale(&style.paneltitle.points, block->aspect);
	ui_fontscale(&style.grouplabel.points, block->aspect);
	ui_fontscale(&style.widgetlabel.points, block->aspect);
	ui_fontscale(&style.widget.points, block->aspect);
	
	/* scale block min/max to rect */
	ui_but_to_pixelrect(&rect, ar, block, NULL);
	
	/* pixel space for AA widgets */
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	
	wmOrtho2(-0.01f, ar->winx - 0.01f, -0.01f, ar->winy - 0.01f);
	
	/* back */
	if (block->flag & UI_BLOCK_LOOP)
		ui_draw_menu_back(&style, block, &rect);
	else if (block->panel)
		ui_draw_aligned_panel(&style, block, &rect);

	/* widgets */
	for (but = block->buttons.first; but; but = but->next) {
		if (!(but->flag & (UI_HIDDEN | UI_SCROLLED))) {
			ui_but_to_pixelrect(&rect, ar, block, but);
		
			/* XXX: figure out why invalid coordinates happen when closing render window */
			/* and material preview is redrawn in main window (temp fix for bug #23848) */
			if (rect.xmin < rect.xmax && rect.ymin < rect.ymax)
				ui_draw_but(C, ar, &style, but, &rect);
		}
	}
	
	/* restore matrix */
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	if (multisample_enabled)
		glEnable(GL_MULTISAMPLE_ARB);
	
	ui_draw_links(block);
}

/* ************* EVENTS ************* */

static void ui_is_but_sel(uiBut *but, double *value)
{
	short is_push = 0;  /* (0 == UNSELECT), (1 == SELECT), (2 == DO-NOHING) */
	short is_true = TRUE;

	if (ELEM3(but->type, TOGN, ICONTOGN, OPTIONN)) {
		is_true = FALSE;
	}

	if (but->bit) {
		int lvalue;
		UI_GET_BUT_VALUE_INIT(but, *value);
		lvalue = (int)*value;
		if (UI_BITBUT_TEST(lvalue, (but->bitnr))) {
			is_push = is_true;
		}
		else {
			is_push = !is_true;
		}
	}
	else {
		switch (but->type) {
			case BUT:
				is_push = 2;
				break;
			case HOTKEYEVT:
			case KEYEVT:
				is_push = 2;
				break;
			case TOGBUT:
			case TOG:
			case TOGR:
			case TOG3:
			case BUT_TOGDUAL:
			case ICONTOG:
			case OPTION:
				UI_GET_BUT_VALUE_INIT(but, *value);
				if (*value != (double)but->hardmin) is_push = 1;
				break;
			case ICONTOGN:
			case TOGN:
			case OPTIONN:
				UI_GET_BUT_VALUE_INIT(but, *value);
				if (*value == 0.0) is_push = 1;
				break;
			case ROW:
			case LISTROW:
				UI_GET_BUT_VALUE_INIT(but, *value);
				/* support for rna enum buts */
				if (but->rnaprop && (RNA_property_flag(but->rnaprop) & PROP_ENUM_FLAG)) {
					if ((int)*value & (int)but->hardmax) is_push = 1;
				}
				else {
					if (*value == (double)but->hardmax) is_push = 1;
				}
				break;
			case COLOR:
				is_push = 2;
				break;
			default:
				is_push = 2;
				break;
		}
	}
	
	if (is_push == 2) ;
	else if (is_push == 1) but->flag |= UI_SELECT;
	else but->flag &= ~UI_SELECT;
}

static uiBut *ui_find_inlink(uiBlock *block, void *poin)
{
	uiBut *but;
	
	but = block->buttons.first;
	while (but) {
		if (but->type == INLINK) {
			if (but->poin == poin) return but;
		}
		but = but->next;
	}
	return NULL;
}

static void ui_add_link_line(ListBase *listb, uiBut *but, uiBut *bt)
{
	uiLinkLine *line;
	
	line = MEM_callocN(sizeof(uiLinkLine), "linkline");
	BLI_addtail(listb, line);
	line->from = but;
	line->to = bt;
}

uiBut *uiFindInlink(uiBlock *block, void *poin)
{
	return ui_find_inlink(block, poin);
}

void uiComposeLinks(uiBlock *block)
{
	uiBut *but, *bt;
	uiLink *link;
	void ***ppoin;
	int a;
	
	but = block->buttons.first;
	while (but) {
		if (but->type == LINK) {
			link = but->link;
			
			/* for all pointers in the array */
			if (link) {
				if (link->ppoin) {
					ppoin = link->ppoin;
					for (a = 0; a < *(link->totlink); a++) {
						bt = ui_find_inlink(block, (*ppoin)[a]);
						if (bt) {
							ui_add_link_line(&link->lines, but, bt);
						}
					}
				}
				else if (link->poin) {
					bt = ui_find_inlink(block, *(link->poin) );
					if (bt) {
						ui_add_link_line(&link->lines, but, bt);
					}
				}
			}
		}
		but = but->next;
	}
}


/* ************************************************ */

void uiBlockSetButLock(uiBlock *block, int val, const char *lockstr)
{
	if (val) {
		block->lock = val ? TRUE : FALSE;
		block->lockstr = lockstr;
	}
}

void uiBlockClearButLock(uiBlock *block)
{
	block->lock = FALSE;
	block->lockstr = NULL;
}

/* *************************************************************** */

void ui_delete_linkline(uiLinkLine *line, uiBut *but)
{
	uiLink *link;
	int a, b;
	
	BLI_remlink(&but->link->lines, line);

	link = line->from->link;

	/* are there more pointers allowed? */
	if (link->ppoin) {
		
		if (*(link->totlink) == 1) {
			*(link->totlink) = 0;
			MEM_freeN(*(link->ppoin));
			*(link->ppoin) = NULL;
		}
		else {
			b = 0;
			for (a = 0; a < (*(link->totlink)); a++) {
				if ((*(link->ppoin))[a] != line->to->poin) {
					(*(link->ppoin))[b] = (*(link->ppoin))[a];
					b++;
				}
			}	
			(*(link->totlink))--;
		}
	}
	else {
		*(link->poin) = NULL;
	}

	MEM_freeN(line);
	//REDRAW
}

/* *********************** data get/set ***********************
 * this either works with the pointed to data, or can work with
 * an edit override pointer while dragging for example */

/* for buttons pointing to color for example */
void ui_get_but_vectorf(uiBut *but, float vec[3])
{
	PropertyRNA *prop;
	int a, tot;

	if (but->editvec) {
		copy_v3_v3(vec, but->editvec);
	}

	if (but->rnaprop) {
		prop = but->rnaprop;

		vec[0] = vec[1] = vec[2] = 0.0f;

		if (RNA_property_type(prop) == PROP_FLOAT) {
			tot = RNA_property_array_length(&but->rnapoin, prop);
			tot = MIN2(tot, 3);

			for (a = 0; a < tot; a++)
				vec[a] = RNA_property_float_get_index(&but->rnapoin, prop, a);
		}
	}
	else if (but->pointype == UI_BUT_POIN_CHAR) {
		char *cp = (char *)but->poin;
		vec[0] = ((float)cp[0]) / 255.0f;
		vec[1] = ((float)cp[1]) / 255.0f;
		vec[2] = ((float)cp[2]) / 255.0f;
	}
	else if (but->pointype == UI_BUT_POIN_FLOAT) {
		float *fp = (float *)but->poin;
		copy_v3_v3(vec, fp);
	}
	else {
		if (but->editvec == NULL) {
			fprintf(stderr, "ui_get_but_vectorf: can't get color, should never happen\n");
			vec[0] = vec[1] = vec[2] = 0.0f;
		}
	}

	if (but->type == BUT_NORMAL) {
		normalize_v3(vec);
	}
}

/* for buttons pointing to color for example */
void ui_set_but_vectorf(uiBut *but, const float vec[3])
{
	PropertyRNA *prop;

	if (but->editvec) {
		copy_v3_v3(but->editvec, vec);
	}

	if (but->rnaprop) {
		prop = but->rnaprop;

		if (RNA_property_type(prop) == PROP_FLOAT) {
			int tot;
			int a;

			tot = RNA_property_array_length(&but->rnapoin, prop);
			tot = MIN2(tot, 3);

			for (a = 0; a < tot; a++) {
				RNA_property_float_set_index(&but->rnapoin, prop, a, vec[a]);
			}
		}
	}
	else if (but->pointype == UI_BUT_POIN_CHAR) {
		char *cp = (char *)but->poin;
		cp[0] = (char)(0.5f + vec[0] * 255.0f);
		cp[1] = (char)(0.5f + vec[1] * 255.0f);
		cp[2] = (char)(0.5f + vec[2] * 255.0f);
	}
	else if (but->pointype == UI_BUT_POIN_FLOAT) {
		float *fp = (float *)but->poin;
		copy_v3_v3(fp, vec);
	}
}

int ui_is_but_float(uiBut *but)
{
	if (but->pointype == UI_BUT_POIN_FLOAT && but->poin)
		return 1;
	
	if (but->rnaprop && RNA_property_type(but->rnaprop) == PROP_FLOAT)
		return 1;
	
	return 0;
}

int ui_is_but_unit(uiBut *but)
{
	UnitSettings *unit = but->block->unit;
	const int unit_type = uiButGetUnitType(but);

	if (unit_type == PROP_UNIT_NONE)
		return 0;

#if 1 /* removed so angle buttons get correct snapping */
	if (unit->system_rotation == USER_UNIT_ROT_RADIANS && unit_type == PROP_UNIT_ROTATION)
		return 0;
#endif
	
	/* for now disable time unit conversion */	
	if (unit_type == PROP_UNIT_TIME)
		return 0;

	if (unit->system == USER_UNIT_NONE) {
		if (unit_type != PROP_UNIT_ROTATION) {
			return 0;
		}
	}

	return 1;
}

int ui_is_but_rna_valid(uiBut *but)
{
	if (but->rnaprop == NULL || RNA_struct_contains_property(&but->rnapoin, but->rnaprop)) {
		return TRUE;
	}
	else {
		printf("property removed %s: %p\n", but->drawstr, but->rnaprop);
		return FALSE;
	}
}

double ui_get_but_val(uiBut *but)
{
	PropertyRNA *prop;
	double value = 0.0;

	if (but->editval) { return *(but->editval); }
	if (but->poin == NULL && but->rnapoin.data == NULL) return 0.0;

	if (but->rnaprop) {
		prop = but->rnaprop;

		switch (RNA_property_type(prop)) {
			case PROP_BOOLEAN:
				if (RNA_property_array_check(prop))
					value = RNA_property_boolean_get_index(&but->rnapoin, prop, but->rnaindex);
				else
					value = RNA_property_boolean_get(&but->rnapoin, prop);
				break;
			case PROP_INT:
				if (RNA_property_array_check(prop))
					value = RNA_property_int_get_index(&but->rnapoin, prop, but->rnaindex);
				else
					value = RNA_property_int_get(&but->rnapoin, prop);
				break;
			case PROP_FLOAT:
				if (RNA_property_array_check(prop))
					value = RNA_property_float_get_index(&but->rnapoin, prop, but->rnaindex);
				else
					value = RNA_property_float_get(&but->rnapoin, prop);
				break;
			case PROP_ENUM:
				value = RNA_property_enum_get(&but->rnapoin, prop);
				break;
			default:
				value = 0.0;
				break;
		}
	}
	else if (but->type == HSVSLI) {
		float *fp, hsv[3];
		
		fp = (but->editvec) ? but->editvec : (float *)but->poin;
		rgb_to_hsv_v(fp, hsv);

		switch (but->str[0]) {
			case 'H': value = hsv[0]; break;
			case 'S': value = hsv[1]; break;
			case 'V': value = hsv[2]; break;
		}
	} 
	else if (but->pointype == UI_BUT_POIN_CHAR) {
		value = *(char *)but->poin;
	}
	else if (but->pointype == UI_BUT_POIN_SHORT) {
		value = *(short *)but->poin;
	} 
	else if (but->pointype == UI_BUT_POIN_INT) {
		value = *(int *)but->poin;
	} 
	else if (but->pointype == UI_BUT_POIN_FLOAT) {
		value = *(float *)but->poin;
	}

	return value;
}

void ui_set_but_val(uiBut *but, double value)
{
	PropertyRNA *prop;

	/* value is a hsv value: convert to rgb */
	if (but->rnaprop) {
		prop = but->rnaprop;

		if (RNA_property_editable(&but->rnapoin, prop)) {
			switch (RNA_property_type(prop)) {
				case PROP_BOOLEAN:
					if (RNA_property_array_length(&but->rnapoin, prop))
						RNA_property_boolean_set_index(&but->rnapoin, prop, but->rnaindex, value);
					else
						RNA_property_boolean_set(&but->rnapoin, prop, value);
					break;
				case PROP_INT:
					if (RNA_property_array_length(&but->rnapoin, prop))
						RNA_property_int_set_index(&but->rnapoin, prop, but->rnaindex, (int)value);
					else
						RNA_property_int_set(&but->rnapoin, prop, (int)value);
					break;
				case PROP_FLOAT:
					if (RNA_property_array_length(&but->rnapoin, prop))
						RNA_property_float_set_index(&but->rnapoin, prop, but->rnaindex, value);
					else
						RNA_property_float_set(&but->rnapoin, prop, value);
					break;
				case PROP_ENUM:
					if (RNA_property_flag(prop) & PROP_ENUM_FLAG) {
						int ivalue = (int)value;
						ivalue ^= RNA_property_enum_get(&but->rnapoin, prop); /* toggle for enum/flag buttons */
						RNA_property_enum_set(&but->rnapoin, prop, ivalue);
					}
					else {
						RNA_property_enum_set(&but->rnapoin, prop, value);
					}
					break;
				default:
					break;
			}
		}

		/* we can't be sure what RNA set functions actually do,
		 * so leave this unset */
		value = UI_BUT_VALUE_UNSET;
	}
	else if (but->pointype == 0) ;
	else if (but->type == HSVSLI) {
		float *fp, hsv[3];
		
		fp = (but->editvec) ? but->editvec : (float *)but->poin;
		rgb_to_hsv_v(fp, hsv);
		
		switch (but->str[0]) {
			case 'H': hsv[0] = value; break;
			case 'S': hsv[1] = value; break;
			case 'V': hsv[2] = value; break;
		}
		
		hsv_to_rgb_v(hsv, fp);
		
	}
	else {
		/* first do rounding */
		if (but->pointype == UI_BUT_POIN_CHAR) {
			value = (char)floor(value + 0.5);
		}
		else if (but->pointype == UI_BUT_POIN_SHORT) {
			/* gcc 3.2.1 seems to have problems
			 * casting a double like 32772.0 to
			 * a short so we cast to an int, then
			 * to a short.
			 *
			 * Update: even in gcc.4.6 using intermediate int cast gives -32764,
			 * where as a direct cast from double to short gives -32768,
			 * if this difference isn't important we could remove this hack,
			 * since we dont support gcc3 anymore - Campbell */
			int gcckludge;
			gcckludge = (int) floor(value + 0.5);
			value = (short)gcckludge;
		}
		else if (but->pointype == UI_BUT_POIN_INT)
			value = (int)floor(value + 0.5);
		else if (but->pointype == UI_BUT_POIN_FLOAT) {
			float fval = (float)value;
			if (fval >= -0.00001f && fval <= 0.00001f) fval = 0.0f;  /* prevent negative zero */
			value = fval;
		}
		
		/* then set value with possible edit override */
		if (but->editval)
			value = *but->editval = value;
		else if (but->pointype == UI_BUT_POIN_CHAR)
			value = *((char *)but->poin) = (char)value;
		else if (but->pointype == UI_BUT_POIN_SHORT)
			value = *((short *)but->poin) = (short)value;
		else if (but->pointype == UI_BUT_POIN_INT)
			value = *((int *)but->poin) = (int)value;
		else if (but->pointype == UI_BUT_POIN_FLOAT)
			value = *((float *)but->poin) = (float)value;
	}

	/* update select flag */
	ui_is_but_sel(but, &value);
}

int ui_get_but_string_max_length(uiBut *but)
{
	if (ELEM(but->type, TEX, SEARCH_MENU))
		return but->hardmax;
	else if (but->type == IDPOIN)
		return MAX_ID_NAME - 2;
	else
		return UI_MAX_DRAW_STR;
}

static double ui_get_but_scale_unit(uiBut *but, double value)
{
	UnitSettings *unit = but->block->unit;
	int unit_type = uiButGetUnitType(but);

	if (unit_type == PROP_UNIT_LENGTH) {
		return value * (double)unit->scale_length;
	}
	else if (unit_type == PROP_UNIT_AREA) {
		return value * pow(unit->scale_length, 2);
	}
	else if (unit_type == PROP_UNIT_VOLUME) {
		return value * pow(unit->scale_length, 3);
	}
	else if (unit_type == PROP_UNIT_TIME) { /* WARNING - using evil_C :| */
		Scene *scene = CTX_data_scene(but->block->evil_C);
		return FRA2TIME(value);
	}
	else {
		return value;
	}
}

/* str will be overwritten */
void ui_convert_to_unit_alt_name(uiBut *but, char *str, size_t maxlen)
{
	if (ui_is_but_unit(but)) {
		UnitSettings *unit = but->block->unit;
		int unit_type = uiButGetUnitType(but);
		char *orig_str;
		
		orig_str = MEM_callocN(sizeof(char) * maxlen + 1, "textedit sub str");
		memcpy(orig_str, str, maxlen);
		
		bUnit_ToUnitAltName(str, maxlen, orig_str, unit->system, RNA_SUBTYPE_UNIT_VALUE(unit_type));
		
		MEM_freeN(orig_str);
	}
}

static void ui_get_but_string_unit(uiBut *but, char *str, int len_max, double value, int pad)
{
	UnitSettings *unit = but->block->unit;
	int do_split = unit->flag & USER_UNIT_OPT_SPLIT;
	int unit_type = uiButGetUnitType(but);
	int precision = but->a2;

	if (unit->scale_length < 0.0001f) unit->scale_length = 1.0f;  // XXX do_versions

	/* Sanity checks */
	if (precision > PRECISION_FLOAT_MAX) precision = PRECISION_FLOAT_MAX;
	else if (precision == 0) precision = 2;

	bUnit_AsString(str, len_max, ui_get_but_scale_unit(but, value), precision,
	               unit->system, RNA_SUBTYPE_UNIT_VALUE(unit_type), do_split, pad);
}

static float ui_get_but_step_unit(uiBut *but, float step_default)
{
	int unit_type = RNA_SUBTYPE_UNIT_VALUE(uiButGetUnitType(but));
	float step;

	step = bUnit_ClosestScalar(ui_get_but_scale_unit(but, step_default), but->block->unit->system, unit_type);

	if (step > 0.0f) { /* -1 is an error value */
		return (float)((double)step / ui_get_but_scale_unit(but, 1.0)) * 100.0f;
	}
	else {
		return step_default;
	}
}


void ui_get_but_string(uiBut *but, char *str, size_t maxlen)
{
	if (but->rnaprop && ELEM3(but->type, TEX, IDPOIN, SEARCH_MENU)) {
		PropertyType type;
		char *buf = NULL;
		int buf_len;

		type = RNA_property_type(but->rnaprop);

		if (type == PROP_STRING) {
			/* RNA string */
			buf = RNA_property_string_get_alloc(&but->rnapoin, but->rnaprop, str, maxlen, &buf_len);
		}
		else if (type == PROP_POINTER) {
			/* RNA pointer */
			PointerRNA ptr = RNA_property_pointer_get(&but->rnapoin, but->rnaprop);
			buf = RNA_struct_name_get_alloc(&ptr, str, maxlen, &buf_len);
		}

		if (!buf) {
			str[0] = '\0';
		}
		else if (buf && buf != str) {
			/* string was too long, we have to truncate */
			memcpy(str, buf, MIN2(maxlen, buf_len + 1));
			MEM_freeN(buf);
		}
	}
	else if (but->type == IDPOIN) {
		/* ID pointer */
		if (but->idpoin_idpp) { /* Can be NULL for ID properties by python */
			ID *id = *(but->idpoin_idpp);
			if (id) {
				BLI_strncpy(str, id->name + 2, maxlen);
				return;
			}
		}
		str[0] = '\0';
		return;
	}
	else if (but->type == TEX) {
		/* string */
		BLI_strncpy(str, but->poin, maxlen);
		return;
	}
	else if (but->type == SEARCH_MENU) {
		/* string */
		BLI_strncpy(str, but->poin, maxlen);
		return;
	}
	else if (ui_but_anim_expression_get(but, str, maxlen))
		;  /* driver expression */
	else {
		/* number editing */
		double value;

		value = ui_get_but_val(but);

		if (ui_is_but_float(but)) {
			if (ui_is_but_unit(but)) {
				ui_get_but_string_unit(but, str, maxlen, value, 0);
			}
			else {
				const int prec = ui_but_float_precision(but, value);
				BLI_snprintf(str, maxlen, "%.*f", prec, value);
			}
		}
		else
			BLI_snprintf(str, maxlen, "%d", (int)value);
	}
}

#ifdef WITH_PYTHON

static int ui_set_but_string_eval_num_unit(bContext *C, uiBut *but, const char *str, double *value)
{
	char str_unit_convert[256];
	const int unit_type = uiButGetUnitType(but);

	BLI_strncpy(str_unit_convert, str, sizeof(str_unit_convert));

	/* ugly, use the draw string to get the value,
	 * this could cause problems if it includes some text which resolves to a unit */
	bUnit_ReplaceString(str_unit_convert, sizeof(str_unit_convert), but->drawstr,
	                    ui_get_but_scale_unit(but, 1.0), but->block->unit->system, RNA_SUBTYPE_UNIT_VALUE(unit_type));

	return (BPY_button_exec(C, str_unit_convert, value, TRUE) != -1);
}

#endif /* WITH_PYTHON */


int ui_set_but_string_eval_num(bContext *C, uiBut *but, const char *str, double *value)
{
	int ok = FALSE;

#ifdef WITH_PYTHON

	if (str[0] != '\0') {
		int is_unit_but = ui_is_but_unit(but);
		/* only enable verbose if we won't run again with units */
		if (BPY_button_exec(C, str, value, is_unit_but == FALSE) != -1) {
			/* if the value parsed ok without unit conversion this button may still need a unit multiplier */
			if (is_unit_but) {
				char str_new[128];

				BLI_snprintf(str_new, sizeof(str_new), "%f", *value);
				ok = ui_set_but_string_eval_num_unit(C, but, str_new, value);
			}
			else {
				ok = TRUE; /* parse normal string via py (no unit conversion needed) */
			}
		}
		else if (is_unit_but) {
			/* parse failed, this is a unit but so run replacements and parse again */
			ok = ui_set_but_string_eval_num_unit(C, but, str, value);
		}
	}

#else /* WITH_PYTHON */

	*value = atof(str);
	ok = TRUE;

	(void)C;
	(void)but;

#endif /* WITH_PYTHON */

	return ok;
}


int ui_set_but_string(bContext *C, uiBut *but, const char *str)
{
	if (but->rnaprop && ELEM3(but->type, TEX, IDPOIN, SEARCH_MENU)) {
		if (RNA_property_editable(&but->rnapoin, but->rnaprop)) {
			PropertyType type;

			type = RNA_property_type(but->rnaprop);

			if (type == PROP_STRING) {
				/* RNA string */
				RNA_property_string_set(&but->rnapoin, but->rnaprop, str);
				return 1;
			}
			else if (type == PROP_POINTER) {
				/* RNA pointer */
				PointerRNA ptr, rptr;
				PropertyRNA *prop;

				if (str == NULL || str[0] == '\0') {
					RNA_property_pointer_set(&but->rnapoin, but->rnaprop, PointerRNA_NULL);
					return 1;
				}
				else {
					ptr = but->rnasearchpoin;
					prop = but->rnasearchprop;
					
					if (prop && RNA_property_collection_lookup_string(&ptr, prop, str, &rptr))
						RNA_property_pointer_set(&but->rnapoin, but->rnaprop, rptr);

					return 1;
				}

				return 0;
			}
		}
	}
	else if (but->type == IDPOIN) {
		/* ID pointer */
		but->idpoin_func(C, str, but->idpoin_idpp);
		return 1;
	}
	else if (but->type == TEX) {
		/* string */
		if (ui_is_but_utf8(but)) BLI_strncpy_utf8(but->poin, str, but->hardmax);
		else BLI_strncpy(but->poin, str, but->hardmax);

		return 1;
	}
	else if (but->type == SEARCH_MENU) {
		/* string */
		BLI_strncpy(but->poin, str, but->hardmax);
		return 1;
	}
	else if (ui_but_anim_expression_set(but, str)) {
		/* driver expression */
		return 1;
	}
	else if (str[0] == '#') {
		/* shortcut to create new driver expression (versus immediate Py-execution) */
		return ui_but_anim_expression_create(but, str + 1);
	}
	else {
		/* number editing */
		double value;

		if (ui_set_but_string_eval_num(C, but, str, &value) == FALSE) {
			return 0;
		}

		if (!ui_is_but_float(but)) value = (int)floor(value + 0.5);
		if (but->type == NUMABS) value = fabs(value);

		/* not that we use hard limits here */
		if (value < (double)but->hardmin) value = but->hardmin;
		if (value > (double)but->hardmax) value = but->hardmax;

		ui_set_but_val(but, value);
		return 1;
	}

	return 0;
}

void ui_set_but_default(bContext *C, short all)
{
	PointerRNA ptr;

	WM_operator_properties_create(&ptr, "UI_OT_reset_default_button");
	RNA_boolean_set(&ptr, "all", all);
	WM_operator_name_call(C, "UI_OT_reset_default_button", WM_OP_EXEC_DEFAULT, &ptr);
	WM_operator_properties_free(&ptr);
}

static double soft_range_round_up(double value, double max)
{
	/* round up to .., 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, .. */
	double newmax = pow(10.0, ceil(log(value) / M_LN10));

	if (newmax * 0.2 >= max && newmax * 0.2 >= value)
		return newmax * 0.2;
	else if (newmax * 0.5 >= max && newmax * 0.5 >= value)
		return newmax * 0.5;
	else
		return newmax;
}

static double soft_range_round_down(double value, double max)
{
	/* round down to .., 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, .. */
	double newmax = pow(10.0, floor(log(value) / M_LN10));

	if (newmax * 5.0 <= max && newmax * 5.0 <= value)
		return newmax * 5.0;
	else if (newmax * 2.0 <= max && newmax * 2.0 <= value)
		return newmax * 2.0;
	else
		return newmax;
}

void ui_set_but_soft_range(uiBut *but, double value)
{
	/* ideally we would not limit this but practically, its more then
	 * enough worst case is very long vectors wont use a smart soft-range
	 * which isn't so bad. */

	if (but->rnaprop) {
		const PropertyType type = RNA_property_type(but->rnaprop);
		double softmin, softmax /*, step, precision*/;
		double value_min = value;
		double value_max = value;

		/* clamp button range to something reasonable in case
		 * we get -inf/inf from RNA properties */
		if (type == PROP_INT) {
			int imin, imax, istep;
			const int array_len = RNA_property_array_length(&but->rnapoin, but->rnaprop);

			RNA_property_int_ui_range(&but->rnapoin, but->rnaprop, &imin, &imax, &istep);
			softmin = (imin == INT_MIN) ? -1e4 : imin;
			softmax = (imin == INT_MAX) ? 1e4 : imax;
			/*step= istep;*/ /*UNUSED*/
			/*precision= 1;*/ /*UNUSED*/

			if (array_len >= 2) {
				int value_range[2];
				RNA_property_int_get_array_range(&but->rnapoin, but->rnaprop, value_range);
				value_min = (double)value_range[0];
				value_max = (double)value_range[1];
			}
		}
		else if (type == PROP_FLOAT) {
			float fmin, fmax, fstep, fprecision;
			const int array_len = RNA_property_array_length(&but->rnapoin, but->rnaprop);

			RNA_property_float_ui_range(&but->rnapoin, but->rnaprop, &fmin, &fmax, &fstep, &fprecision);
			softmin = (fmin == -FLT_MAX) ? (float)-1e4 : fmin;
			softmax = (fmax == FLT_MAX) ? (float)1e4 : fmax;
			/*step= fstep;*/ /*UNUSED*/
			/*precision= fprecision;*/ /*UNUSED*/

			if (array_len >= 2) {
				float value_range[2];
				RNA_property_float_get_array_range(&but->rnapoin, but->rnaprop, value_range);
				value_min = (double)value_range[0];
				value_max = (double)value_range[1];
			}
		}
		else
			return;

		/* if the value goes out of the soft/max range, adapt the range */
		if (value_min + 1e-10 < softmin) {
			if (value_min < 0.0)
				softmin = -soft_range_round_up(-value_min, -softmin);
			else
				softmin = soft_range_round_down(value_min, softmin);

			if (softmin < (double)but->hardmin)
				softmin = (double)but->hardmin;
		}
		if (value_max - 1e-10 > softmax) {
			if (value_max < 0.0)
				softmax = -soft_range_round_down(-value_max, -softmax);
			else
				softmax = soft_range_round_up(value_max, softmax);

			if (softmax > (double)but->hardmax)
				softmax = but->hardmax;
		}

		but->softmin = softmin;
		but->softmax = softmax;
	}
}

/* ******************* Free ********************/

static void ui_free_link(uiLink *link)
{
	if (link) {	
		BLI_freelistN(&link->lines);
		MEM_freeN(link);
	}
}

/* can be called with C==NULL */
static void ui_free_but(const bContext *C, uiBut *but)
{
	if (but->opptr) {
		WM_operator_properties_free(but->opptr);
		MEM_freeN(but->opptr);
	}

	if (but->func_argN) {
		MEM_freeN(but->func_argN);
	}

	if (but->active) {
		/* XXX solve later, buttons should be free-able without context ideally,
		 * however they may have open tooltips or popup windows, which need to
		 * be closed using a context pointer */
		if (C) {
			ui_button_active_free(C, but);
		}
		else {
			if (but->active) {
				MEM_freeN(but->active);
			}
		}
	}
	if (but->str && but->str != but->strdata) {
		MEM_freeN(but->str);
	}
	ui_free_link(but->link);

	if ((but->type == BUT_IMAGE) && but->poin) {
		IMB_freeImBuf((struct ImBuf *)but->poin);
	}

	MEM_freeN(but);
}

/* can be called with C==NULL */
void uiFreeBlock(const bContext *C, uiBlock *block)
{
	uiBut *but;

	while ( (but = block->buttons.first) ) {
		BLI_remlink(&block->buttons, but);	
		ui_free_but(C, but);
	}

	if (block->unit) {
		MEM_freeN(block->unit);
	}

	if (block->func_argN) {
		MEM_freeN(block->func_argN);
	}

	CTX_store_free_list(&block->contexts);

	BLI_freelistN(&block->saferct);
	
	MEM_freeN(block);
}

/* can be called with C==NULL */
void uiFreeBlocks(const bContext *C, ListBase *lb)
{
	uiBlock *block;
	
	while ( (block = lb->first) ) {
		BLI_remlink(lb, block);
		uiFreeBlock(C, block);
	}
}

void uiFreeInactiveBlocks(const bContext *C, ListBase *lb)
{
	uiBlock *block, *nextblock;

	for (block = lb->first; block; block = nextblock) {
		nextblock = block->next;
	
		if (!block->handle) {
			if (!block->active) {
				BLI_remlink(lb, block);
				uiFreeBlock(C, block);
			}
			else
				block->active = 0;
		}
	}
}

void uiBlockSetRegion(uiBlock *block, ARegion *region)
{
	ListBase *lb = &region->uiblocks;
	uiBlock *oldblock = NULL;

	/* each listbase only has one block with this name, free block
	 * if is already there so it can be rebuilt from scratch */
	if (lb) {
		oldblock = BLI_findstring(lb, block->name, offsetof(uiBlock, name));

		if (oldblock) {
			oldblock->active = 0;
			oldblock->panel = NULL;
		}

		/* at the beginning of the list! for dynamical menus/blocks */
		BLI_addhead(lb, block);
	}

	block->oldblock = oldblock;
}

uiBlock *uiBeginBlock(const bContext *C, ARegion *region, const char *name, short dt)
{
	uiBlock *block;
	wmWindow *window;
	Scene *scn;
	int getsizex, getsizey;

	window = CTX_wm_window(C);
	scn = CTX_data_scene(C);

	block = MEM_callocN(sizeof(uiBlock), "uiBlock");
	block->active = 1;
	block->dt = dt;
	block->evil_C = (void *)C;  /* XXX */

	if (scn) {
		block->color_profile = TRUE;

		/* store display device name, don't lookup for transformations yet
		 * block could be used for non-color displays where looking up for transformation
		 * would slow down redraw, so only lookup for actual transform when it's indeed
		 * needed
		 */
		block->display_device = scn->display_settings.display_device;

		/* copy to avoid crash when scene gets deleted with ui still open */
		block->unit = MEM_mallocN(sizeof(scn->unit), "UI UnitSettings");
		memcpy(block->unit, &scn->unit, sizeof(scn->unit));
	}

	BLI_strncpy(block->name, name, sizeof(block->name));

	if (region)
		uiBlockSetRegion(block, region);

	/* window matrix and aspect */
	if (region && region->swinid) {
		wm_subwindow_getmatrix(window, region->swinid, block->winmat);
		wm_subwindow_getsize(window, region->swinid, &getsizex, &getsizey);

		/* TODO - investigate why block->winmat[0][0] is negative
		 * in the image view when viewRedrawForce is called */
		block->aspect = 2.0f / fabsf(getsizex * block->winmat[0][0]);
	}
	else {
		/* no subwindow created yet, for menus for example, so we
		 * use the main window instead, since buttons are created
		 * there anyway */
		wm_subwindow_getmatrix(window, window->screen->mainwin, block->winmat);
		wm_subwindow_getsize(window, window->screen->mainwin, &getsizex, &getsizey);

		block->aspect = 2.0f / fabsf(getsizex * block->winmat[0][0]);
		block->auto_open = TRUE;
		block->flag |= UI_BLOCK_LOOP; /* tag as menu */
	}

	return block;
}

uiBlock *uiGetBlock(const char *name, ARegion *ar)
{
	return BLI_findstring(&ar->uiblocks, name, offsetof(uiBlock, name));
}

void uiBlockSetEmboss(uiBlock *block, char dt)
{
	block->dt = dt;
}

void ui_check_but(uiBut *but)
{
	/* if something changed in the button */
	double value = UI_BUT_VALUE_UNSET;
//	float okwidth; // UNUSED
	
	ui_is_but_sel(but, &value);
	
	/* only update soft range while not editing */
	if (but->rnaprop && !(but->editval || but->editstr || but->editvec)) {
		UI_GET_BUT_VALUE_INIT(but, value);
		ui_set_but_soft_range(but, value);
	}

	/* test for min and max, icon sliders, etc */
	switch (but->type) {
		case NUM:
		case SLI:
		case SCROLL:
		case NUMSLI:
		case HSVSLI:
			UI_GET_BUT_VALUE_INIT(but, value);
			if (value < (double)but->hardmin) ui_set_but_val(but, but->hardmin);
			else if (value > (double)but->hardmax) ui_set_but_val(but, but->hardmax);
			break;
			
		case NUMABS:
		{
			double value_abs;
			UI_GET_BUT_VALUE_INIT(but, value);
			value_abs = fabs(value);
			if (value_abs < (double)but->hardmin) ui_set_but_val(but, but->hardmin);
			else if (value_abs > (double)but->hardmax) ui_set_but_val(but, but->hardmax);
			break;
		}
		case ICONTOG: 
		case ICONTOGN:
			if (!but->rnaprop || (RNA_property_flag(but->rnaprop) & PROP_ICONS_CONSECUTIVE)) {
				if (but->flag & UI_SELECT) but->iconadd = 1;
				else but->iconadd = 0;
			}
			break;
			
		case ICONROW:
			if (!but->rnaprop || (RNA_property_flag(but->rnaprop) & PROP_ICONS_CONSECUTIVE)) {
				UI_GET_BUT_VALUE_INIT(but, value);
				but->iconadd = (int)value - (int)(but->hardmin);
			}
			break;
			
		case ICONTEXTROW:
			if (!but->rnaprop || (RNA_property_flag(but->rnaprop) & PROP_ICONS_CONSECUTIVE)) {
				UI_GET_BUT_VALUE_INIT(but, value);
				but->iconadd = (int)value - (int)(but->hardmin);
			}
			break;

			/* quiet warnings for unhandled types */
		default:
			break;
	}
	
	
	/* safety is 4 to enable small number buttons (like 'users') */
	// okwidth= -4 + (BLI_rcti_size_x(&but->rect)); // UNUSED
	
	/* name: */
	switch (but->type) {
	
		case MENU:
		case ICONTEXTROW:
		
			if (BLI_rctf_size_x(&but->rect) > 24.0f) {
				UI_GET_BUT_VALUE_INIT(but, value);
				ui_set_name_menu(but, (int)value);
			}
			break;
	
		case NUM:
		case NUMSLI:
		case HSVSLI:
		case NUMABS:

			UI_GET_BUT_VALUE_INIT(but, value);

			if (ui_is_but_float(but)) {
				if (value == (double) FLT_MAX) {
					BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%sinf", but->str);
				}
				else if (value == (double) -FLT_MAX) {
					BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%s-inf", but->str);
				}
				/* support length type buttons */
				else if (ui_is_but_unit(but)) {
					char new_str[sizeof(but->drawstr)];
					ui_get_but_string_unit(but, new_str, sizeof(new_str), value, TRUE);
					BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%s%s", but->str, new_str);
				}
				else {
					const int prec = ui_but_float_precision(but, value);
					BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%s%.*f", but->str, prec, value);
				}
			}
			else {
				BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%s%d", but->str, (int)value);
			}
			
			if (but->rnaprop) {
				PropertySubType pstype = RNA_property_subtype(but->rnaprop);
			
				if (pstype == PROP_PERCENTAGE)
					strcat(but->drawstr, "%");
			}
			break;

		case LABEL:
			if (ui_is_but_float(but)) {
				int prec;
				UI_GET_BUT_VALUE_INIT(but, value);
				prec = ui_but_float_precision(but, value);
				BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%s%.*f", but->str, prec, value);
			}
			else {
				BLI_strncpy(but->drawstr, but->str, UI_MAX_DRAW_STR);
			}
		
			break;

		case IDPOIN:
		case TEX:
		case SEARCH_MENU:
			if (!but->editstr) {
				char str[UI_MAX_DRAW_STR];

				ui_get_but_string(but, str, UI_MAX_DRAW_STR - strlen(but->str));

				BLI_snprintf(but->drawstr, sizeof(but->drawstr), "%s%s", but->str, str);
			}
			break;
	
		case KEYEVT:
			BLI_strncpy(but->drawstr, but->str, UI_MAX_DRAW_STR);
			if (but->flag & UI_SELECT) {
				strcat(but->drawstr, "Press a key");
			}
			else {
				UI_GET_BUT_VALUE_INIT(but, value);
				strcat(but->drawstr, WM_key_event_string((short)value));
			}
			break;
		
		case HOTKEYEVT:
			if (but->flag & UI_SELECT) {
				but->drawstr[0] = '\0';

				if (but->modifier_key) {
					char *str = but->drawstr;

					if (but->modifier_key & KM_SHIFT)
						str = strcat(str, "Shift ");
					if (but->modifier_key & KM_CTRL)
						str = strcat(str, "Ctrl ");
					if (but->modifier_key & KM_ALT)
						str = strcat(str, "Alt ");
					if (but->modifier_key & KM_OSKEY)
						str = strcat(str, "Cmd ");

					(void)str; /* UNUSED */
				}
				else
					strcat(but->drawstr, "Press a key  ");
			}
			else
				BLI_strncpy(but->drawstr, but->str, UI_MAX_DRAW_STR);

			break;

		case BUT_TOGDUAL:
			/* trying to get the dual-icon to left of text... not very nice */
			if (but->str[0]) {
				BLI_strncpy(but->drawstr, "  ", UI_MAX_DRAW_STR);
				BLI_strncpy(but->drawstr + 2, but->str, UI_MAX_DRAW_STR - 2);
			}
			break;
		
		case HSVCUBE:
		case HSVCIRCLE:
			break;
		default:
			BLI_strncpy(but->drawstr, but->str, UI_MAX_DRAW_STR);
		
	}

	/* if we are doing text editing, this will override the drawstr */
	if (but->editstr)
		BLI_strncpy(but->drawstr, but->editstr, UI_MAX_DRAW_STR);
	
	/* text clipping moved to widget drawing code itself */
}


void uiBlockBeginAlign(uiBlock *block)
{
	/* if other align was active, end it */
	if (block->flag & UI_BUT_ALIGN) uiBlockEndAlign(block);

	block->flag |= UI_BUT_ALIGN_DOWN;	
	block->alignnr++;

	/* buttons declared after this call will get this align nr */ // XXX flag?
}

static int buts_are_horiz(uiBut *but1, uiBut *but2)
{
	float dx, dy;
	
	dx = fabs(but1->rect.xmax - but2->rect.xmin);
	dy = fabs(but1->rect.ymin - but2->rect.ymax);
	
	if (dx > dy) return 0;
	return 1;
}

void uiBlockEndAlign(uiBlock *block)
{
	block->flag &= ~UI_BUT_ALIGN;  /* all 4 flags */
}

int ui_but_can_align(uiBut *but)
{
	return !ELEM3(but->type, LABEL, OPTION, OPTIONN);
}

static void ui_block_do_align_but(uiBut *first, short nr)
{
	uiBut *prev, *but = NULL, *next;
	int flag = 0, cols = 0, rows = 0;
	
	/* auto align */

	for (but = first; but && but->alignnr == nr; but = but->next) {
		if (but->next && but->next->alignnr == nr) {
			if (buts_are_horiz(but, but->next)) cols++;
			else rows++;
		}
	}

	/* rows==0: 1 row, cols==0: 1 column */
	
	/* note;  how it uses 'flag' in loop below (either set it, or OR it) is confusing */
	for (but = first, prev = NULL; but && but->alignnr == nr; prev = but, but = but->next) {
		next = but->next;
		if (next && next->alignnr != nr)
			next = NULL;

		/* clear old flag */
		but->flag &= ~UI_BUT_ALIGN;
			
		if (flag == 0) {  /* first case */
			if (next) {
				if (buts_are_horiz(but, next)) {
					if (rows == 0)
						flag = UI_BUT_ALIGN_RIGHT;
					else 
						flag = UI_BUT_ALIGN_DOWN | UI_BUT_ALIGN_RIGHT;
				}
				else {
					flag = UI_BUT_ALIGN_DOWN;
				}
			}
		}
		else if (next == NULL) {  /* last case */
			if (prev) {
				if (buts_are_horiz(prev, but)) {
					if (rows == 0)
						flag = UI_BUT_ALIGN_LEFT;
					else
						flag = UI_BUT_ALIGN_TOP | UI_BUT_ALIGN_LEFT;
				}
				else flag = UI_BUT_ALIGN_TOP;
			}
		}
		else if (buts_are_horiz(but, next)) {
			/* check if this is already second row */
			if (prev && buts_are_horiz(prev, but) == 0) {
				flag &= ~UI_BUT_ALIGN_LEFT;
				flag |= UI_BUT_ALIGN_TOP;
				/* exception case: bottom row */
				if (rows > 0) {
					uiBut *bt = but;
					while (bt && bt->alignnr == nr) {
						if (bt->next && bt->next->alignnr == nr && buts_are_horiz(bt, bt->next) == 0) break;
						bt = bt->next;
					}
					if (bt == NULL || bt->alignnr != nr) flag = UI_BUT_ALIGN_TOP | UI_BUT_ALIGN_RIGHT;
				}
			}
			else flag |= UI_BUT_ALIGN_LEFT;
		}
		else {
			if (cols == 0) {
				flag |= UI_BUT_ALIGN_TOP;
			}
			else {  /* next button switches to new row */
				
				if (prev && buts_are_horiz(prev, but))
					flag |= UI_BUT_ALIGN_LEFT;
				else {
					flag &= ~UI_BUT_ALIGN_LEFT;
					flag |= UI_BUT_ALIGN_TOP;
				}
				
				if ((flag & UI_BUT_ALIGN_TOP) == 0) {  /* stil top row */
					if (prev) {
						if (next && buts_are_horiz(but, next))
							flag = UI_BUT_ALIGN_DOWN | UI_BUT_ALIGN_LEFT | UI_BUT_ALIGN_RIGHT;
						else {
							/* last button in top row */
							flag = UI_BUT_ALIGN_DOWN | UI_BUT_ALIGN_LEFT;
						}
					}
					else 
						flag |= UI_BUT_ALIGN_DOWN;
				}
				else 
					flag |= UI_BUT_ALIGN_TOP;
			}
		}
		
		but->flag |= flag;
		
		/* merge coordinates */
		if (prev) {
			/* simple cases */
			if (rows == 0) {
				but->rect.xmin = (prev->rect.xmax + but->rect.xmin) / 2.0f;
				prev->rect.xmax = but->rect.xmin;
			}
			else if (cols == 0) {
				but->rect.ymax = (prev->rect.ymin + but->rect.ymax) / 2.0f;
				prev->rect.ymin = but->rect.ymax;
			}
			else {
				if (buts_are_horiz(prev, but)) {
					but->rect.xmin = (prev->rect.xmax + but->rect.xmin) / 2.0f;
					prev->rect.xmax = but->rect.xmin;
					/* copy height too */
					but->rect.ymax = prev->rect.ymax;
				}
				else if (prev->prev && buts_are_horiz(prev->prev, prev) == 0) {
					/* the previous button is a single one in its row */
					but->rect.ymax = (prev->rect.ymin + but->rect.ymax) / 2.0f;
					prev->rect.ymin = but->rect.ymax;
					
					but->rect.xmin = prev->rect.xmin;
					if (next && buts_are_horiz(but, next) == 0)
						but->rect.xmax = prev->rect.xmax;
				}
				else {
					/* the previous button is not a single one in its row */
					but->rect.ymax = prev->rect.ymin;
				}
			}
		}
	}
}

void ui_block_do_align(uiBlock *block)
{
	uiBut *but;
	short nr;

	/* align buttons with same align nr */
	for (but = block->buttons.first; but; ) {
		if (but->alignnr) {
			nr = but->alignnr;
			ui_block_do_align_but(but, nr);

			/* skip with same number */
			for (; but && but->alignnr == nr; but = but->next) {
				/* pass */
			}

			if (!but) {
				break;
			}
		}
		else {
			but = but->next;
		}
	}
}

struct ColorManagedDisplay *ui_block_display_get(uiBlock *block)
{
	return IMB_colormanagement_display_get_named(block->display_device);
}

void ui_block_to_display_space_v3(uiBlock *block, float pixel[3])
{
	struct ColorManagedDisplay *display = ui_block_display_get(block);

	IMB_colormanagement_scene_linear_to_display_v3(pixel, display);
}

void ui_block_to_scene_linear_v3(uiBlock *block, float pixel[3])
{
	struct ColorManagedDisplay *display = ui_block_display_get(block);

	IMB_colormanagement_display_to_scene_linear_v3(pixel, display);
}

/**
 * \brief ui_def_but is the function that draws many button types
 *
 * \param x,y The lower left hand corner of the button (X axis)
 * \param width,height The size of the button.
 *
 * for float buttons:
 * - \a a1 Click Step (how much to change the value each click)
 * - \a a2 Number of decimal point values to display. 0 defaults to 3 (0.000)
 *      1,2,3, and a maximum of 4, all greater values will be clamped to 4.
 */
static uiBut *ui_def_but(uiBlock *block, int type, int retval, const char *str,
                         int x, int y, short width, short height,
                         void *poin, float min, float max, float a1, float a2, const char *tip)
{
	uiBut *but;
	int slen;

	/* we could do some more error checks here */
	if ((type & BUTTYPE) == LABEL) {
		BLI_assert((poin != NULL || min != 0.0f || max != 0.0f || (a1 == 0.0f && a2 != 0.0f) || (a1 != 0.0f && a1 != 1.0f)) == FALSE);
	}

	if (type & UI_BUT_POIN_TYPES) {  /* a pointer is required */
		if (poin == NULL) {
			BLI_assert(0);
			return NULL;
		}
	}

	but = MEM_callocN(sizeof(uiBut), "uiBut");

	but->type = type & BUTTYPE;
	but->pointype = type & UI_BUT_POIN_TYPES;
	but->bit = type & UI_BUT_POIN_BIT;
	but->bitnr = type & 31;
	but->icon = ICON_NONE;
	but->iconadd = 0;

	but->retval = retval;

	slen = strlen(str);
	if (slen >= UI_MAX_NAME_STR - 1) {
		but->str = MEM_mallocN(slen + 2, "ui_def_but str"); /* why +2 ? */
	}
	else {
		but->str = but->strdata;
	}
	memcpy(but->str, str, slen + 1);

	but->rect.xmin = x;
	but->rect.ymin = y;
	but->rect.xmax = but->rect.xmin + width;
	but->rect.ymax = but->rect.ymin + height;

	but->poin = poin;
	but->hardmin = but->softmin = min;
	but->hardmax = but->softmax = max;
	but->a1 = a1;
	but->a2 = a2;
	but->tip = tip;

	but->lock = block->lock;
	but->lockstr = block->lockstr;
	but->dt = block->dt;

	but->aspect = 1.0f;  /* XXX block->aspect; */
	but->block = block;  /* pointer back, used for frontbuffer status, and picker */

	if ((block->flag & UI_BUT_ALIGN) && ui_but_can_align(but))
		but->alignnr = block->alignnr;

	but->func = block->func;
	but->func_arg1 = block->func_arg1;
	but->func_arg2 = block->func_arg2;

	but->funcN = block->funcN;
	if (block->func_argN)
		but->func_argN = MEM_dupallocN(block->func_argN);
	
	but->pos = -1;   /* cursor invisible */

	if (ELEM4(but->type, NUM, NUMABS, NUMSLI, HSVSLI)) {    /* add a space to name */
		/* slen remains unchanged from previous assignment, ensure this stays true */
		if (slen > 0 && slen < UI_MAX_NAME_STR - 2) {
			if (but->str[slen - 1] != ' ') {
				but->str[slen] = ' ';
				but->str[slen + 1] = 0;
			}
		}
	}

	if ((block->flag & UI_BLOCK_LOOP) ||
	    ELEM8(but->type, MENU, TEX, LABEL, IDPOIN, BLOCK, BUTM, SEARCH_MENU, PROGRESSBAR))
	{
		but->flag |= (UI_TEXT_LEFT | UI_ICON_LEFT);
	}
	else if (but->type == BUT_TOGDUAL) {
		but->flag |= UI_ICON_LEFT;
	}

	but->flag |= (block->flag & UI_BUT_ALIGN);

	if (but->lock == TRUE) {
		if (but->lockstr) {
			but->flag |= UI_BUT_DISABLED;
		}
	}

	/* keep track of UI_interface.h */
	if (ELEM7(but->type, BLOCK, BUT, LABEL, PULLDOWN, ROUNDBOX, LISTBOX, BUTM)) ;
	else if (ELEM(but->type, SCROLL, SEPR /* , FTPREVIEW */ )) ;
	else if (but->type >= SEARCH_MENU) ;
	else but->flag |= UI_BUT_UNDO;

	BLI_addtail(&block->buttons, but);
	
	if (block->curlayout)
		ui_layout_add_but(block->curlayout, but);

#ifdef WITH_PYTHON
	/* if the 'UI_OT_editsource' is running, extract the source info from the button  */
	if (UI_editsource_enable_check()) {
		UI_editsource_active_but_test(but);
	}
#endif

	return but;
}

/* ui_def_but_rna_propname and ui_def_but_rna
 * both take the same args except for propname vs prop, this is done so we can
 * avoid an extra lookup on 'prop' when its already available.
 *
 * When this kind of change won't disrupt branches, best look into making more
 * of our UI functions take prop rather then propname.
 */

#define UI_DEF_BUT_RNA_DISABLE(but)  { \
		but->flag |= UI_BUT_DISABLED;  \
		but->lock = TRUE;              \
		but->lockstr = "";             \
	} (void)0


static uiBut *ui_def_but_rna(uiBlock *block, int type, int retval, const char *str,
							 int x, int y, short width, short height,
                             PointerRNA *ptr, PropertyRNA *prop, int index,
                             float min, float max, float a1, float a2,  const char *tip)
{
	const PropertyType proptype = RNA_property_type(prop);
	uiBut *but;
	int freestr = 0, icon = 0;

	/* use rna values if parameters are not specified */
	if (!str) {
		if (type == MENU && proptype == PROP_ENUM) {
			EnumPropertyItem *item;
			DynStr *dynstr;
			int i, totitem, value, free;

			RNA_property_enum_items_gettexted(block->evil_C, ptr, prop, &item, &totitem, &free);
			value = RNA_property_enum_get(ptr, prop);

			dynstr = BLI_dynstr_new();
			BLI_dynstr_appendf(dynstr, "%s%%t", RNA_property_ui_name(prop));
			for (i = 0; i < totitem; i++) {
				if (!item[i].identifier[0]) {
					if (item[i].name)
						BLI_dynstr_appendf(dynstr, "|%s%%l", item[i].name);
					else
						BLI_dynstr_append(dynstr, "|%l");
				}
				else if (item[i].icon)
					BLI_dynstr_appendf(dynstr, "|%s %%i%d %%x%d", item[i].name, item[i].icon, item[i].value);
				else
					BLI_dynstr_appendf(dynstr, "|%s %%x%d", item[i].name, item[i].value);

				if (value == item[i].value)
					icon = item[i].icon;
			}
			str = BLI_dynstr_get_cstring(dynstr);
			BLI_dynstr_free(dynstr);

			if (free) {
				MEM_freeN(item);
			}

			freestr = 1;
		}
		else if (ELEM(type, ROW, LISTROW) && proptype == PROP_ENUM) {
			EnumPropertyItem *item;
			int i, totitem, free;

			RNA_property_enum_items_gettexted(block->evil_C, ptr, prop, &item, &totitem, &free);
			for (i = 0; i < totitem; i++) {
				if (item[i].identifier[0] && item[i].value == (int)max) {
					str = item[i].name;
					icon = item[i].icon;
				}
			}

			if (!str) {
				str = RNA_property_ui_name(prop);
			}
			if (free) {
				MEM_freeN(item);
			}
		}
		else {
			str = RNA_property_ui_name(prop);
			icon = RNA_property_ui_icon(prop);
		}
	}

	if (!tip && proptype != PROP_ENUM)
		tip = RNA_property_ui_description(prop);

	if (min == max || a1 == -1 || a2 == -1) {
		if (proptype == PROP_INT) {
			int hardmin, hardmax, softmin, softmax, step;

			RNA_property_int_range(ptr, prop, &hardmin, &hardmax);
			RNA_property_int_ui_range(ptr, prop, &softmin, &softmax, &step);

			if (!ELEM(type, ROW, LISTROW) && min == max) {
				min = hardmin;
				max = hardmax;
			}
			if (a1 == -1)
				a1 = step;
			if (a2 == -1)
				a2 = 0;
		}
		else if (proptype == PROP_FLOAT) {
			float hardmin, hardmax, softmin, softmax, step, precision;

			RNA_property_float_range(ptr, prop, &hardmin, &hardmax);
			RNA_property_float_ui_range(ptr, prop, &softmin, &softmax, &step, &precision);

			if (!ELEM(type, ROW, LISTROW) && min == max) {
				min = hardmin;
				max = hardmax;
			}
			if (a1 == -1)
				a1 = step;
			if (a2 == -1)
				a2 = precision;
		}
		else if (proptype == PROP_STRING) {
			min = 0;
			max = RNA_property_string_maxlength(prop);
			if (max == 0) /* interface code should ideally support unlimited length */
				max = UI_MAX_DRAW_STR;
		}
	}

	/* now create button */
	but = ui_def_but(block, type, retval, str, x, y, width, height, NULL, min, max, a1, a2, tip);

	but->rnapoin = *ptr;
	but->rnaprop = prop;

	if (RNA_property_array_length(&but->rnapoin, but->rnaprop))
		but->rnaindex = index;
	else
		but->rnaindex = 0;

	if (icon) {
		but->icon = (BIFIconID)icon;
		but->flag |= UI_HAS_ICON;
		but->flag |= UI_ICON_LEFT;
	}
	
	if (!RNA_property_editable(&but->rnapoin, prop)) {
		UI_DEF_BUT_RNA_DISABLE(but);
	}

	if (but->flag & UI_BUT_UNDO && (ui_but_is_rna_undo(but) == FALSE)) {
		but->flag &= ~UI_BUT_UNDO;
	}

	/* If this button uses units, calculate the step from this */
	if ((proptype == PROP_FLOAT) && ui_is_but_unit(but)) {
		but->a1 = ui_get_but_step_unit(but, but->a1);
	}

	if (freestr) {
		MEM_freeN((void *)str);
	}

	return but;
}

static uiBut *ui_def_but_rna_propname(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, PointerRNA *ptr, const char *propname, int index, float min, float max, float a1, float a2,  const char *tip)
{
	PropertyRNA *prop = RNA_struct_find_property(ptr, propname);
	uiBut *but;

	if (prop) {
		but = ui_def_but_rna(block, type, retval, str, x, y, width, height, ptr, prop, index, min, max, a1, a2,  tip);
	}
	else {
		but = ui_def_but(block, type, retval, propname, x, y, width, height, NULL, min, max, a1, a2, tip);

		UI_DEF_BUT_RNA_DISABLE(but);
	}

	return but;
}

static uiBut *ui_def_but_operator_ptr(uiBlock *block, int type, wmOperatorType *ot, int opcontext, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but;

	if (!str) {
		if (ot && ot->srna)
			str = RNA_struct_ui_name(ot->srna);
		else
			str = "";
	}

	if ((!tip || tip[0] == '\0') && ot && ot->srna) {
		tip = RNA_struct_ui_description(ot->srna);
	}

	but = ui_def_but(block, type, -1, str, x, y, width, height, NULL, 0, 0, 0, 0, tip);
	but->optype = ot;
	but->opcontext = opcontext;
	but->flag &= ~UI_BUT_UNDO; /* no need for ui_but_is_undo(), we never need undo here */

	if (!ot) {
		but->flag |= UI_BUT_DISABLED;
		but->lock = TRUE;
		but->lockstr = "";
	}

	return but;
}

#if 0 /* UNUSED */
static uiBut *UNUSED_FUNCTION(ui_def_but_operator) (uiBlock * block, int type, const char *opname, int opcontext, const char *str, int x, int y, short width, short height, const char *tip)
{
	wmOperatorType *ot = WM_operatortype_find(opname, 0);
	if (str == NULL && ot == NULL) str = opname;
	return ui_def_but_operator_ptr(block, type, ot, opcontext, str, x, y, width, height, tip);
}
#endif

static uiBut *ui_def_but_operator_text(uiBlock *block, int type, const char *opname, int opcontext, const char *str, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2, const char *tip)
{
	uiBut *but;
	wmOperatorType *ot;
	
	ot = WM_operatortype_find(opname, 0);

	if (!str) {
		if (ot) str = ot->name;
		else str = opname;
	}
	
	if ((!tip || tip[0] == '\0') && ot && ot->description) {
		tip = ot->description;
	}

	but = ui_def_but(block, type, -1, str, x, y, width, height, poin, min, max, a1, a2, tip);
	but->optype = ot;
	but->opcontext = opcontext;
	but->flag &= ~UI_BUT_UNDO; /* no need for ui_but_is_undo(), we never need undo here */

	if (!ot) {
		but->flag |= UI_BUT_DISABLED;
		but->lock = TRUE;
		but->lockstr = "";
	}

	return but;
}

uiBut *uiDefBut(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2, const char *tip)
{
	uiBut *but = ui_def_but(block, type, retval, str, x, y, width, height, poin, min, max, a1, a2, tip);

	ui_check_but(but);
	
	return but;
}

/* if _x_ is a power of two (only one bit) return the power,
 * otherwise return -1.
 * (1<<findBitIndex(x))==x for powers of two.
 */
static int findBitIndex(unsigned int x)
{
	if (!x || !is_power_of_2_i(x)) { /* is_power_of_2_i(x) strips lowest bit */
		return -1;
	}
	else {
		int idx = 0;

		if (x & 0xFFFF0000) idx += 16, x >>= 16;
		if (x & 0xFF00) idx += 8, x >>= 8;
		if (x & 0xF0) idx += 4, x >>= 4;
		if (x & 0xC) idx += 2, x >>= 2;
		if (x & 0x2) idx += 1;

		return idx;
	}
}

/* autocomplete helper functions */
struct AutoComplete {
	size_t maxlen;
	char *truncate;
	const char *startname;
};

AutoComplete *autocomplete_begin(const char *startname, size_t maxlen)
{
	AutoComplete *autocpl;
	
	autocpl = MEM_callocN(sizeof(AutoComplete), "AutoComplete");
	autocpl->maxlen = maxlen;
	autocpl->truncate = MEM_callocN(sizeof(char) * maxlen, "AutoCompleteTruncate");
	autocpl->startname = startname;

	return autocpl;
}

void autocomplete_do_name(AutoComplete *autocpl, const char *name)
{
	char *truncate = autocpl->truncate;
	const char *startname = autocpl->startname;
	int a;

	for (a = 0; a < autocpl->maxlen - 1; a++) {
		if (startname[a] == 0 || startname[a] != name[a])
			break;
	}
	/* found a match */
	if (startname[a] == 0) {
		/* first match */
		if (truncate[0] == 0)
			BLI_strncpy(truncate, name, autocpl->maxlen);
		else {
			/* remove from truncate what is not in bone->name */
			for (a = 0; a < autocpl->maxlen - 1; a++) {
				if (name[a] == 0) {
					truncate[a] = 0;
					break;
				}
				else if (truncate[a] != name[a])
					truncate[a] = 0;
			}
		}
	}
}

void autocomplete_end(AutoComplete *autocpl, char *autoname)
{	
	if (autocpl->truncate[0])
		BLI_strncpy(autoname, autocpl->truncate, autocpl->maxlen);
	else {
		if (autoname != autocpl->startname) /* don't copy a string over its self */
			BLI_strncpy(autoname, autocpl->startname, autocpl->maxlen);
	}
	MEM_freeN(autocpl->truncate);
	MEM_freeN(autocpl);
}

/* autocomplete callback for ID buttons */
static void autocomplete_id(bContext *C, char *str, void *arg_v)
{
	int blocktype = (intptr_t)arg_v;
	ListBase *listb = which_libbase(CTX_data_main(C), blocktype);
	
	if (listb == NULL) return;
	
	/* search if str matches the beginning of an ID struct */
	if (str[0]) {
		AutoComplete *autocpl = autocomplete_begin(str, MAX_ID_NAME - 2);
		ID *id;
		
		for (id = listb->first; id; id = id->next)
			autocomplete_do_name(autocpl, id->name + 2);

		autocomplete_end(autocpl, str);
	}
}

static void ui_check_but_and_iconize(uiBut *but, int icon)
{
	if (icon) {
		but->icon = (BIFIconID) icon;
		but->flag |= UI_HAS_ICON;
	}

	ui_check_but(but);
}

static uiBut *uiDefButBit(uiBlock *block, int type, int bit, int retval, const char *str, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2,  const char *tip)
{
	int bitIdx = findBitIndex(bit);
	if (bitIdx == -1) {
		return NULL;
	}
	else {
		return uiDefBut(block, type | UI_BUT_POIN_BIT | bitIdx, retval, str, x, y, width, height, poin, min, max, a1, a2, tip);
	}
}
uiBut *uiDefButF(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, float *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefBut(block, type | UI_BUT_POIN_FLOAT, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButBitF(uiBlock *block, int type, int bit, int retval, const char *str, int x, int y, short width, short height, float *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefButBit(block, type | UI_BUT_POIN_FLOAT, bit, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButI(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, int *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefBut(block, type | UI_BUT_POIN_INT, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButBitI(uiBlock *block, int type, int bit, int retval, const char *str, int x, int y, short width, short height, int *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefButBit(block, type | UI_BUT_POIN_INT, bit, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButS(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, short *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefBut(block, type | UI_BUT_POIN_SHORT, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButBitS(uiBlock *block, int type, int bit, int retval, const char *str, int x, int y, short width, short height, short *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefButBit(block, type | UI_BUT_POIN_SHORT, bit, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButC(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, char *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefBut(block, type | UI_BUT_POIN_CHAR, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButBitC(uiBlock *block, int type, int bit, int retval, const char *str, int x, int y, short width, short height, char *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefButBit(block, type | UI_BUT_POIN_CHAR, bit, retval, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefButR(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, PointerRNA *ptr, const char *propname, int index, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but;
	but = ui_def_but_rna_propname(block, type, retval, str, x, y, width, height, ptr, propname, index, min, max, a1, a2, tip);
	ui_check_but(but);
	return but;
}
uiBut *uiDefButR_prop(uiBlock *block, int type, int retval, const char *str, int x, int y, short width, short height, PointerRNA *ptr, PropertyRNA *prop, int index, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but;
	but = ui_def_but_rna(block, type, retval, str, x, y, width, height, ptr, prop, index, min, max, a1, a2, tip);
	ui_check_but(but);
	return but;
}

uiBut *uiDefButO_ptr(uiBlock *block, int type, wmOperatorType *ot, int opcontext, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but;
	but = ui_def_but_operator_ptr(block, type, ot, opcontext, str, x, y, width, height, tip);
	ui_check_but(but);
	return but;
}
uiBut *uiDefButO(uiBlock *block, int type, const char *opname, int opcontext, const char *str, int x, int y, short width, short height, const char *tip)
{
	wmOperatorType *ot = WM_operatortype_find(opname, 0);
	if (str == NULL && ot == NULL) str = opname;
	return uiDefButO_ptr(block, type, ot, opcontext, str, x, y, width, height, tip);
}

uiBut *uiDefButTextO(uiBlock *block, int type, const char *opname, int opcontext, const char *str, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but = ui_def_but_operator_text(block, type, opname, opcontext, str, x, y, width, height, poin, min, max, a1, a2, tip);
	ui_check_but(but);
	return but;
}

/* if a1==1.0 then a2 is an extra icon blending factor (alpha 0.0 - 1.0) */
uiBut *uiDefIconBut(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but = ui_def_but(block, type, retval, "", x, y, width, height, poin, min, max, a1, a2, tip);
	ui_check_but_and_iconize(but, icon);
	return but;
}
static uiBut *uiDefIconButBit(uiBlock *block, int type, int bit, int retval, int icon, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2, const char *tip)
{
	int bitIdx = findBitIndex(bit);
	if (bitIdx == -1) {
		return NULL;
	}
	else {
		return uiDefIconBut(block, type | UI_BUT_POIN_BIT | bitIdx, retval, icon, x, y, width, height, poin, min, max, a1, a2, tip);
	}
}

uiBut *uiDefIconButF(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, float *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconBut(block, type | UI_BUT_POIN_FLOAT, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButBitF(uiBlock *block, int type, int bit, int retval, int icon, int x, int y, short width, short height, float *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconButBit(block, type | UI_BUT_POIN_FLOAT, bit, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButI(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, int *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconBut(block, type | UI_BUT_POIN_INT, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButBitI(uiBlock *block, int type, int bit, int retval, int icon, int x, int y, short width, short height, int *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconButBit(block, type | UI_BUT_POIN_INT, bit, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButS(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, short *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconBut(block, type | UI_BUT_POIN_SHORT, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButBitS(uiBlock *block, int type, int bit, int retval, int icon, int x, int y, short width, short height, short *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconButBit(block, type | UI_BUT_POIN_SHORT, bit, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButC(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, char *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconBut(block, type | UI_BUT_POIN_CHAR, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButBitC(uiBlock *block, int type, int bit, int retval, int icon, int x, int y, short width, short height, char *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconButBit(block, type | UI_BUT_POIN_CHAR, bit, retval, icon, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconButR(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, PointerRNA *ptr, const char *propname, int index, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but;
	but = ui_def_but_rna_propname(block, type, retval, "", x, y, width, height, ptr, propname, index, min, max, a1, a2, tip);
	ui_check_but_and_iconize(but, icon);
	return but;
}
uiBut *uiDefIconButR_prop(uiBlock *block, int type, int retval, int icon, int x, int y, short width, short height, PointerRNA *ptr, PropertyRNA *prop, int index, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but;
	but = ui_def_but_rna(block, type, retval, "", x, y, width, height, ptr, prop, index, min, max, a1, a2, tip);
	ui_check_but_and_iconize(but, icon);
	return but;
}

uiBut *uiDefIconButO_ptr(uiBlock *block, int type, wmOperatorType *ot, int opcontext, int icon, int x, int y, short width, short height, const char *tip)
{
	uiBut *but;
	but = ui_def_but_operator_ptr(block, type, ot, opcontext, "", x, y, width, height, tip);
	ui_check_but_and_iconize(but, icon);
	return but;
}
uiBut *uiDefIconButO(uiBlock *block, int type, const char *opname, int opcontext, int icon, int x, int y, short width, short height, const char *tip)
{
	wmOperatorType *ot = WM_operatortype_find(opname, 0);
	return uiDefIconButO_ptr(block, type, ot, opcontext, icon, x, y, width, height, tip);
}

/* Button containing both string label and icon */
uiBut *uiDefIconTextBut(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but = ui_def_but(block, type, retval, str, x, y, width, height, poin, min, max, a1, a2, tip);
	ui_check_but_and_iconize(but, icon);
	but->flag |= UI_ICON_LEFT;
	return but;
}
static uiBut *uiDefIconTextButBit(uiBlock *block, int type, int bit, int retval, int icon, const char *str, int x, int y, short width, short height, void *poin, float min, float max, float a1, float a2,  const char *tip)
{
	int bitIdx = findBitIndex(bit);
	if (bitIdx == -1) {
		return NULL;
	}
	else {
		return uiDefIconTextBut(block, type | UI_BUT_POIN_BIT | bitIdx, retval, icon, str, x, y, width, height, poin, min, max, a1, a2, tip);
	}
}

uiBut *uiDefIconTextButF(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, float *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextBut(block, type | UI_BUT_POIN_FLOAT, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButBitF(uiBlock *block, int type, int bit, int retval, int icon, const char *str, int x, int y, short width, short height, float *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextButBit(block, type | UI_BUT_POIN_FLOAT, bit, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButI(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, int *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextBut(block, type | UI_BUT_POIN_INT, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButBitI(uiBlock *block, int type, int bit, int retval, int icon, const char *str, int x, int y, short width, short height, int *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextButBit(block, type | UI_BUT_POIN_INT, bit, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButS(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, short *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextBut(block, type | UI_BUT_POIN_SHORT, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButBitS(uiBlock *block, int type, int bit, int retval, int icon, const char *str, int x, int y, short width, short height, short *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextButBit(block, type | UI_BUT_POIN_SHORT, bit, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButC(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, char *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextBut(block, type | UI_BUT_POIN_CHAR, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButBitC(uiBlock *block, int type, int bit, int retval, int icon, const char *str, int x, int y, short width, short height, char *poin, float min, float max, float a1, float a2,  const char *tip)
{
	return uiDefIconTextButBit(block, type | UI_BUT_POIN_CHAR, bit, retval, icon, str, x, y, width, height, (void *) poin, min, max, a1, a2, tip);
}
uiBut *uiDefIconTextButR(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, PointerRNA *ptr, const char *propname, int index, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but;
	but = ui_def_but_rna_propname(block, type, retval, str, x, y, width, height, ptr, propname, index, min, max, a1, a2, tip);
	ui_check_but_and_iconize(but, icon);
	but->flag |= UI_ICON_LEFT;
	return but;
}
uiBut *uiDefIconTextButR_prop(uiBlock *block, int type, int retval, int icon, const char *str, int x, int y, short width, short height, PointerRNA *ptr, PropertyRNA *prop, int index, float min, float max, float a1, float a2,  const char *tip)
{
	uiBut *but;
	but = ui_def_but_rna(block, type, retval, str, x, y, width, height, ptr, prop, index, min, max, a1, a2, tip);
	ui_check_but_and_iconize(but, icon);
	but->flag |= UI_ICON_LEFT;
	return but;
}
uiBut *uiDefIconTextButO_ptr(uiBlock *block, int type, wmOperatorType *ot, int opcontext, int icon, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but;
	but = ui_def_but_operator_ptr(block, type, ot, opcontext, str, x, y, width, height, tip);
	ui_check_but_and_iconize(but, icon);
	but->flag |= UI_ICON_LEFT;
	return but;
}
uiBut *uiDefIconTextButO(uiBlock *block, int type, const char *opname, int opcontext, int icon, const char *str, int x, int y, short width, short height, const char *tip)
{
	wmOperatorType *ot = WM_operatortype_find(opname, 0);
	return uiDefIconTextButO_ptr(block, type, ot, opcontext, icon, str, x, y, width, height, tip);
}

/* END Button containing both string label and icon */

void uiSetButLink(uiBut *but, void **poin, void ***ppoin, short *tot, int from, int to)
{
	uiLink *link;
	
	link = but->link = MEM_callocN(sizeof(uiLink), "new uilink");
	
	link->poin = poin;
	link->ppoin = ppoin;
	link->totlink = tot;
	link->fromcode = from;
	link->tocode = to;
}

/* cruft to make uiBlock and uiBut private */

int uiBlocksGetYMin(ListBase *lb)
{
	uiBlock *block;
	int min = 0;
	
	for (block = lb->first; block; block = block->next)
		if (block == lb->first || block->rect.ymin < min)
			min = block->rect.ymin;
			
	return min;
}

void uiBlockSetDirection(uiBlock *block, int direction)
{
	block->direction = direction;
}

/* this call escapes if there's alignment flags */
void uiBlockFlipOrder(uiBlock *block)
{
	ListBase lb;
	uiBut *but, *next;
	float centy, miny = 10000, maxy = -10000;

	if (U.uiflag & USER_MENUFIXEDORDER)
		return;
	else if (block->flag & UI_BLOCK_NO_FLIP)
		return;
	
	for (but = block->buttons.first; but; but = but->next) {
		if (but->flag & UI_BUT_ALIGN) return;
		if (but->rect.ymin < miny) miny = but->rect.ymin;
		if (but->rect.ymax > maxy) maxy = but->rect.ymax;
	}
	/* mirror trick */
	centy = (miny + maxy) / 2.0f;
	for (but = block->buttons.first; but; but = but->next) {
		but->rect.ymin = centy - (but->rect.ymin - centy);
		but->rect.ymax = centy - (but->rect.ymax - centy);
		SWAP(float, but->rect.ymin, but->rect.ymax);
	}
	
	/* also flip order in block itself, for example for arrowkey */
	lb.first = lb.last = NULL;
	but = block->buttons.first;
	while (but) {
		next = but->next;
		BLI_remlink(&block->buttons, but);
		BLI_addtail(&lb, but);
		but = next;
	}
	block->buttons = lb;
}


void uiBlockSetFlag(uiBlock *block, int flag)
{
	block->flag |= flag;
}

void uiBlockClearFlag(uiBlock *block, int flag)
{
	block->flag &= ~flag;
}

void uiBlockSetXOfs(uiBlock *block, int xofs)
{
	block->xofs = xofs;
}

void uiButSetFlag(uiBut *but, int flag)
{
	but->flag |= flag;
}

void uiButClearFlag(uiBut *but, int flag)
{
	but->flag &= ~flag;
}

void uiButSetDrawFlag(uiBut *but, int flag)
{
	but->drawflag |= flag;
}

void uiButClearDrawFlag(uiBut *but, int flag)
{
	but->drawflag &= ~flag;
}

int uiButGetRetVal(uiBut *but)
{
	return but->retval;
}

void uiButSetDragID(uiBut *but, ID *id)
{
	but->dragtype = WM_DRAG_ID;
	but->dragpoin = (void *)id;
}

void uiButSetDragRNA(uiBut *but, PointerRNA *ptr)
{
	but->dragtype = WM_DRAG_RNA;
	but->dragpoin = (void *)ptr;
}

void uiButSetDragPath(uiBut *but, const char *path)
{
	but->dragtype = WM_DRAG_PATH;
	but->dragpoin = (void *)path;
}

void uiButSetDragName(uiBut *but, const char *name)
{
	but->dragtype = WM_DRAG_NAME;
	but->dragpoin = (void *)name;
}

/* value from button itself */
void uiButSetDragValue(uiBut *but)
{
	but->dragtype = WM_DRAG_VALUE;
}

void uiButSetDragImage(uiBut *but, const char *path, int icon, struct ImBuf *imb, float scale)
{
	but->dragtype = WM_DRAG_PATH;
	but->icon = icon; /* no flag UI_HAS_ICON, so icon doesnt draw in button */
	but->dragpoin = (void *)path;
	but->imb = imb;
	but->imb_scale = scale;
}

PointerRNA *uiButGetOperatorPtrRNA(uiBut *but)
{
	if (but->optype && !but->opptr) {
		but->opptr = MEM_callocN(sizeof(PointerRNA), "uiButOpPtr");
		WM_operator_properties_create_ptr(but->opptr, but->optype);
	}

	return but->opptr;
}

void uiButSetUnitType(uiBut *but, const int unit_type)
{
	but->unit_type = (unsigned char)(RNA_SUBTYPE_UNIT_VALUE(unit_type));
}

int uiButGetUnitType(uiBut *but)
{
	int ownUnit = (int)but->unit_type;
	
	/* own unit define always takes precedence over RNA provided, allowing for overriding 
	 * default value provided in RNA in a few special cases (i.e. Active Keyframe in Graph Edit)
	 */
	/* XXX: this doesn't allow clearing unit completely, though the same could be said for icons */
	if ((ownUnit != 0) || (but->rnaprop == NULL)) {
		return ownUnit << 16;
	}
	else {
		return RNA_SUBTYPE_UNIT(RNA_property_subtype(but->rnaprop));
	}
}

void uiBlockSetHandleFunc(uiBlock *block, uiBlockHandleFunc func, void *arg)
{
	block->handle_func = func;
	block->handle_func_arg = arg;
}

void uiBlockSetButmFunc(uiBlock *block, uiMenuHandleFunc func, void *arg)
{
	block->butm_func = func;
	block->butm_func_arg = arg;
}

void uiBlockSetFunc(uiBlock *block, uiButHandleFunc func, void *arg1, void *arg2)
{
	block->func = func;
	block->func_arg1 = arg1;
	block->func_arg2 = arg2;
}

void uiBlockSetNFunc(uiBlock *block, uiButHandleFunc func, void *argN, void *arg2)
{
	if (block->func_argN) {
		MEM_freeN(block->func_argN);
	}

	block->funcN = func;
	block->func_argN = argN;
	block->func_arg2 = arg2;
}

void uiButSetRenameFunc(uiBut *but, uiButHandleRenameFunc func, void *arg1)
{
	but->rename_func = func;
	but->rename_arg1 = arg1;
}

void uiBlockSetDrawExtraFunc(uiBlock *block, void (*func)(const bContext *C, void *idv, void *arg1, void *arg2, rcti *rect), void *arg1, void *arg2)
{
	block->drawextra = func;
	block->drawextra_arg1 = arg1;
	block->drawextra_arg2 = arg2;
}

void uiButSetFunc(uiBut *but, uiButHandleFunc func, void *arg1, void *arg2)
{
	but->func = func;
	but->func_arg1 = arg1;
	but->func_arg2 = arg2;
}

void uiButSetNFunc(uiBut *but, uiButHandleNFunc funcN, void *argN, void *arg2)
{
	if (but->func_argN) {
		MEM_freeN(but->func_argN);
	}

	but->funcN = funcN;
	but->func_argN = argN;
	but->func_arg2 = arg2;
}

void uiButSetCompleteFunc(uiBut *but, uiButCompleteFunc func, void *arg)
{
	but->autocomplete_func = func;
	but->autofunc_arg = arg;
}

uiBut *uiDefIDPoinBut(uiBlock *block, uiIDPoinFuncFP func, short blocktype, int retval, const char *str, int x, int y, short width, short height, void *idpp, const char *tip)
{
	uiBut *but = ui_def_but(block, IDPOIN, retval, str, x, y, width, height, NULL, 0.0, 0.0, 0.0, 0.0, tip);
	but->idpoin_func = func;
	but->idpoin_idpp = (ID **) idpp;
	ui_check_but(but);
	
	if (blocktype)
		uiButSetCompleteFunc(but, autocomplete_id, (void *)(intptr_t)blocktype);

	return but;
}

uiBut *uiDefBlockBut(uiBlock *block, uiBlockCreateFunc func, void *arg, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, BLOCK, 0, str, x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);
	but->block_create_func = func;
	ui_check_but(but);
	return but;
}

uiBut *uiDefBlockButN(uiBlock *block, uiBlockCreateFunc func, void *argN, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, BLOCK, 0, str, x, y, width, height, NULL, 0.0, 0.0, 0.0, 0.0, tip);
	but->block_create_func = func;
	if (but->func_argN) {
		MEM_freeN(but->func_argN);
	}
	but->func_argN = argN;
	ui_check_but(but);
	return but;
}


uiBut *uiDefPulldownBut(uiBlock *block, uiBlockCreateFunc func, void *arg, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, PULLDOWN, 0, str, x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);
	but->block_create_func = func;
	ui_check_but(but);
	return but;
}

uiBut *uiDefMenuBut(uiBlock *block, uiMenuCreateFunc func, void *arg, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, PULLDOWN, 0, str, x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);
	but->menu_create_func = func;
	ui_check_but(but);
	return but;
}

uiBut *uiDefIconTextMenuBut(uiBlock *block, uiMenuCreateFunc func, void *arg, int icon, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, PULLDOWN, 0, str, x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);

	but->icon = (BIFIconID) icon;
	but->flag |= UI_HAS_ICON;

	but->flag |= UI_ICON_LEFT;
	but->flag |= UI_ICON_SUBMENU;

	but->menu_create_func = func;
	ui_check_but(but);

	return but;
}

uiBut *uiDefIconMenuBut(uiBlock *block, uiMenuCreateFunc func, void *arg, int icon, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, PULLDOWN, 0, "", x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);

	but->icon = (BIFIconID) icon;
	but->flag |= UI_HAS_ICON;
	but->flag &= ~UI_ICON_LEFT;

	but->menu_create_func = func;
	ui_check_but(but);

	return but;
}

/* Block button containing both string label and icon */
uiBut *uiDefIconTextBlockBut(uiBlock *block, uiBlockCreateFunc func, void *arg, int icon, const char *str, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, BLOCK, 0, str, x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);
	
	/* XXX temp, old menu calls pass on icon arrow, which is now UI_ICON_SUBMENU flag */
	if (icon != ICON_RIGHTARROW_THIN) {
		but->icon = (BIFIconID) icon;
		but->flag |= UI_ICON_LEFT;
	}
	but->flag |= UI_HAS_ICON;
	but->flag |= UI_ICON_SUBMENU;

	but->block_create_func = func;
	ui_check_but(but);
	
	return but;
}

/* Block button containing icon */
uiBut *uiDefIconBlockBut(uiBlock *block, uiBlockCreateFunc func, void *arg, int retval, int icon, int x, int y, short width, short height, const char *tip)
{
	uiBut *but = ui_def_but(block, BLOCK, retval, "", x, y, width, height, arg, 0.0, 0.0, 0.0, 0.0, tip);
	
	but->icon = (BIFIconID) icon;
	but->flag |= UI_HAS_ICON;
	
	but->flag |= UI_ICON_LEFT;
	
	but->block_create_func = func;
	ui_check_but(but);
	
	return but;
}

uiBut *uiDefKeyevtButS(uiBlock *block, int retval, const char *str, int x, int y, short width, short height, short *spoin, const char *tip)
{
	uiBut *but = ui_def_but(block, KEYEVT | UI_BUT_POIN_SHORT, retval, str, x, y, width, height, spoin, 0.0, 0.0, 0.0, 0.0, tip);
	ui_check_but(but);
	return but;
}

/* short pointers hardcoded */
/* modkeypoin will be set to KM_SHIFT, KM_ALT, KM_CTRL, KM_OSKEY bits */
uiBut *uiDefHotKeyevtButS(uiBlock *block, int retval, const char *str, int x, int y, short width, short height, short *keypoin, short *modkeypoin, const char *tip)
{
	uiBut *but = ui_def_but(block, HOTKEYEVT | UI_BUT_POIN_SHORT, retval, str, x, y, width, height, keypoin, 0.0, 0.0, 0.0, 0.0, tip);
	but->modifier_key = *modkeypoin;
	ui_check_but(but);
	return but;
}


/* arg is pointer to string/name, use uiButSetSearchFunc() below to make this work */
/* here a1 and a2, if set, control thumbnail preview rows/cols */
uiBut *uiDefSearchBut(uiBlock *block, void *arg, int retval, int icon, int maxlen, int x, int y, short width, short height, float a1, float a2, const char *tip)
{
	uiBut *but = ui_def_but(block, SEARCH_MENU, retval, "", x, y, width, height, arg, 0.0, maxlen, a1, a2, tip);
	
	but->icon = (BIFIconID) icon;
	but->flag |= UI_HAS_ICON;
	
	but->flag |= UI_ICON_LEFT | UI_TEXT_LEFT;
	
	ui_check_but(but);
	
	return but;
}


/* arg is user value, searchfunc and handlefunc both get it as arg */
/* if active set, button opens with this item visible and selected */
void uiButSetSearchFunc(uiBut *but, uiButSearchFunc sfunc, void *arg, uiButHandleFunc bfunc, void *active)
{
	but->search_func = sfunc;
	but->search_arg = arg;
	
	uiButSetFunc(but, bfunc, arg, active);
	
	/* search buttons show red-alert if item doesn't exist, not for menus */
	if (0 == (but->block->flag & UI_BLOCK_LOOP)) {
		/* skip empty buttons, not all buttons need input, we only show invalid */
		if (but->drawstr[0])
			ui_but_search_test(but);
	}
}

/* push a new event onto event queue to activate the given button 
 * (usually a text-field) upon entering a popup
 */
void uiButSetFocusOnEnter(wmWindow *win, uiBut *but)
{
	wmEvent event;
	
	event = *(win->eventstate);
	event.type = EVT_BUT_OPEN;
	event.val = KM_PRESS;
	event.customdata = but;
	event.customdatafree = FALSE;
	
	wm_event_add(win, &event);
}

void uiButGetStrInfo(bContext *C, uiBut *but, int nbr, ...)
{
	va_list args;

	EnumPropertyItem *items = NULL, *item = NULL;
	int totitems, free_items = FALSE;

	va_start(args, nbr);
	while (nbr--) {
		uiStringInfo *si = (uiStringInfo *) va_arg(args, void *);
		int type = si->type;
		char *tmp = NULL;

		if (type == BUT_GET_LABEL) {
			if (but->str) {
				/* Menu labels can have some complex formating stuff marked by pipes or %t, we don't want those here! */
				const char *tc1, *tc2;

				tc1 = strstr(but->str, "%t");
				tc2 = strstr(but->str, "|"); /* XXX For some reason strchr seems to not work here? */

				if (tc2 && (!tc1 || tc1 > tc2))
					tc1 = tc2;

				if (tc1)
					tmp = BLI_strdupn(but->str, tc1 - but->str);
				else
					tmp = BLI_strdup(but->str);
			}
			else
				type = BUT_GET_RNA_LABEL;  /* Fail-safe solution... */
		}
		else if (type == BUT_GET_TIP) {
			if (but->tip && but->tip[0])
				tmp = BLI_strdup(but->tip);
			else
				type = BUT_GET_RNA_TIP;  /* Fail-safe solution... */
		}

		if (type == BUT_GET_RNAPROP_IDENTIFIER) {
			if (but->rnaprop)
				tmp = BLI_strdup(RNA_property_identifier(but->rnaprop));
		}
		else if (type == BUT_GET_RNASTRUCT_IDENTIFIER) {
			if (but->rnaprop)
				tmp = BLI_strdup(RNA_struct_identifier(but->rnapoin.type));
			else if (but->optype)
				tmp = BLI_strdup(but->optype->idname);
			else if (ELEM(but->type, MENU, PULLDOWN)) {
				MenuType *mt = uiButGetMenuType(but);
				if (mt)
					tmp = BLI_strdup(mt->idname);
			}
		}
		else if (ELEM(type, BUT_GET_RNA_LABEL, BUT_GET_RNA_TIP)) {
			if (but->rnaprop) {
				if (type == BUT_GET_RNA_LABEL)
					tmp = BLI_strdup(RNA_property_ui_name(but->rnaprop));
				else {
					const char *t = RNA_property_ui_description(but->rnaprop);
					if (t && t[0])
						tmp = BLI_strdup(t);
				}
			}
			else if (but->optype) {
				if (type == BUT_GET_RNA_LABEL)
					tmp = BLI_strdup(RNA_struct_ui_name(but->optype->srna));
				else {
					const char *t = RNA_struct_ui_description(but->optype->srna);
					if (t && t[0])
						tmp = BLI_strdup(t);
				}
			}
			else if (ELEM(but->type, MENU, PULLDOWN)) {
				MenuType *mt = uiButGetMenuType(but);
				if (mt) {
					/* not all menus are from python */
					if (mt->ext.srna) {
						if (type == BUT_GET_RNA_LABEL)
							tmp = BLI_strdup(RNA_struct_ui_name(mt->ext.srna));
						else {
							const char *t = RNA_struct_ui_description(mt->ext.srna);
							if (t && t[0])
								tmp = BLI_strdup(t);
						}
					}
				}
			}
		}
		else if (type == BUT_GET_RNA_LABEL_CONTEXT) {
			if (but->rnaprop)
				tmp = BLI_strdup(RNA_property_translation_context(but->rnaprop));
			else if (but->optype)
				tmp = BLI_strdup(RNA_struct_translation_context(but->optype->srna));
			else if (ELEM(but->type, MENU, PULLDOWN)) {
				MenuType *mt = uiButGetMenuType(but);
				if (mt)
					tmp = BLI_strdup(RNA_struct_translation_context(mt->ext.srna));
			}
		}
		else if (ELEM3(type, BUT_GET_RNAENUM_IDENTIFIER, BUT_GET_RNAENUM_LABEL, BUT_GET_RNAENUM_TIP)) {
			if (but->rnaprop && RNA_property_type(but->rnaprop) == PROP_ENUM) {
				if (!item) {
					int i;
					int value = (but->type == ROW) ? (int)but->hardmax : (int)ui_get_but_val(but);
					RNA_property_enum_items_gettexted(C, &but->rnapoin, but->rnaprop, &items, &totitems, &free_items);

					for (i = 0, item = items; i < totitems; i++, item++) {
						if (item->identifier[0] && item->value == value)
							break;
					}
				}
				if (item && item->identifier) {
					if (type == BUT_GET_RNAENUM_IDENTIFIER)
						tmp = BLI_strdup(item->identifier);
					else if (type == BUT_GET_RNAENUM_LABEL)
						tmp = BLI_strdup(item->name);
					else if (item->description && item->description[0])
						tmp = BLI_strdup(item->description);
				}
			}
		}
		else if (type == BUT_GET_OP_KEYMAP) {
			if (but->optype && !(but->block->flag & UI_BLOCK_LOOP)) {
				/* operator keymap (not menus, they already have it) */
				IDProperty *prop = (but->opptr) ? but->opptr->data : NULL;
				char buf[512];

				if (WM_key_event_operator_string(C, but->optype->idname, but->opcontext, prop, TRUE,
				                                 buf, sizeof(buf)))
				{
					tmp = BLI_strdup(buf);
				}
			}
		}

		si->strinfo = tmp;
	}

	if (free_items && items)
		MEM_freeN(items);
}

/* Program Init/Exit */

void UI_init(void)
{
	ui_resources_init();
}

/* after reading userdef file */
void UI_init_userdef(void)
{
	/* fix saved themes */
	init_userdef_do_versions();
	uiStyleInit();
}

void UI_reinit_font(void)
{
	uiStyleInit();
}

void UI_exit(void)
{
	ui_resources_free();
}


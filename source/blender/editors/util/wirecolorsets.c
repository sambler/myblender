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
 * The Original Code is Copyright (C) 2009 Blender Foundation, Shane Ambler
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Shane Ambler
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/util/wirecolorsets.c
 *  \ingroup edutil
 */

 
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_math.h"
#include "BLI_dynstr.h"
#include "BLI_utildefines.h"

#include "DNA_anim_types.h"
#include "DNA_scene_types.h"
#include "DNA_object_types.h"

#include "BKE_main.h"
#include "BKE_animsys.h"
#include "BKE_context.h"
#include "BKE_depsgraph.h"
#include "BKE_report.h"

#include "ED_keyframing.h"
#include "ED_screen.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "WM_api.h"
#include "WM_types.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

/* ************************************************** */
/* WIRECOLOR SETS - OPERATORS (for use in UI panels) */
/* These operators are really duplication of existing functionality, but just for completeness,
 * they're here too, and will give the basic data needed...
 */

/* poll callback for adding default WirecolorSet */
static int wirecolorset_poll_default_add (bContext *C)
{
	/* as long as there's an active Scene, it's fine */
	return (CTX_data_scene(C) != NULL);
}

/* poll callback for editing active WirecolorSet */
static int wirecolorset_poll_active_edit (bContext *C)
{
	Scene *scene= CTX_data_scene(C);
	
	if (scene == NULL)
		return 0;
	
	/* there must be an active WirecolorSet (and WirecolorSets) */
	return ((scene->active_wirecolorset > 0) && (scene->wirecolorsets.first));
}
 
/* Add a Default (Empty) Wirecolor Set ------------------------- */

static int add_default_wirecolorset_exec (bContext *C, wmOperator *UNUSED(op))
{
	Scene *scene= CTX_data_scene(C);
	short flag=0, keyingflag=0;
	
	/* call the API func, and set the active keyingset index */
	BKE_wirecolorset_add(&scene->wirecolorsets, NULL);
	
	scene->active_wirecolorset= BLI_countlist(&scene->wirecolorsets);
	
	/* send notifiers */
	WM_event_add_notifier(C, NC_SCENE|ND_WIRECOLORSET, NULL);
	
	return OPERATOR_FINISHED;
}

void SCENE_OT_wirecolor_set_add (wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Add Empty Wirecolor Set";
	ot->idname= "SCENE_OT_wirecolor_set_add";
	ot->description= "Add a new (empty) Wirecolor Set to the active Scene";
	
	/* callbacks */
	ot->exec= add_default_wirecolorset_exec;
	ot->poll= wirecolorset_poll_default_add;
}

/* Remove 'Active' Wirecolor Set ------------------------- */

static int remove_active_wirecolorset_exec (bContext *C, wmOperator *op)
{
	Scene *scene= CTX_data_scene(C);
	WirecolorSet *wcs;
	
	/* verify the Wirecolor Set to use:
	 *	- use the active one
	 *	- return error if it doesn't exist
	 */
	if (scene->active_wirecolorset == 0) {
		BKE_report(op->reports, RPT_ERROR, "No active Wirecolor Set to remove");
		return OPERATOR_CANCELLED;
	}
	else
		wcs= BLI_findlink(&scene->wirecolorsets, scene->active_wirecolorset-1);
	
	/* free WirecolorSet's data, then remove it from the scene */
	BKE_keyingset_free(wcs);
	BLI_freelinkN(&scene->wirecolorsets, wcs);
	
	/* the active one should now be the previously second-to-last one */
	scene->active_wirecolorset--;
	
	/* send notifiers */
	WM_event_add_notifier(C, NC_SCENE|ND_WIRECOLORSET, NULL);
	
	return OPERATOR_FINISHED;
}

void SCENE_OT_wirecolor_set_remove (wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Removed Active Wirecolor Set";
	ot->idname= "SCENE_OT_wirecolor_set_remove";
	ot->description= "Remove the active Wirecolor Set";
	
	/* callbacks */
	ot->exec= remove_active_wirecolorset_exec;
	ot->poll= wirecolorset_poll_active_edit;
}

/* ******************************************* */

/* Change Active WirecolorSet Operator ------------------------ */
/* This operator checks if a menu should be shown for choosing the WirecolorSet to make the active one */

static int wirecolorset_active_menu_invoke (bContext *C, wmOperator *op, wmEvent *UNUSED(event))
{
	/* call the menu, which will call this operator again, hence the cancelled */
	SCENE_wirecolor_sets_menu_setup(C, op->type->name, "SCENE_OT_wirecolor_set_active_set");
	return OPERATOR_CANCELLED;
}

static int wirecolorset_active_menu_exec (bContext *C, wmOperator *op)
{
	Scene *scene= CTX_data_scene(C);
	int type= RNA_int_get(op->ptr, "type");
	
	/* simply set the scene's active wirecolor set index, unless the type == 0
	 * (i.e. which happens if we want the current active to be maintained) 
	 */
	if (type)
		scene->active_wirecolorset= type;
		
	/* send notifiers */
	WM_event_add_notifier(C, NC_SCENE|ND_WIRECOLORSET, NULL);
	
	return OPERATOR_FINISHED;
}
 
void SCENE_OT_wirecolor_set_active_set (wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Set Active Wirecolor Set";
	ot->idname= "SCENE_OT_wirecolor_set_active_set";
	
	/* callbacks */
	ot->invoke= wirecolorset_active_menu_invoke;
	ot->exec= wirecolorset_active_menu_exec;
	ot->poll= ED_operator_areaactive;
	
	/* flags */
	ot->flag= OPTYPE_REGISTER|OPTYPE_UNDO;
	
	/* wirecolorset to use
	 *	- here the type is int not enum, since many of the indices here are determined dynamically
	 */
	RNA_def_int(ot->srna, "type", 0, INT_MIN, INT_MAX, "Wirecolor Set Number", "Index (determined internally) of the Wirecolor Set to use", 0, 1);
}

/* ******************************************* */
/* WIRECOLOR SETS API (for UI) */

/* Getters for Active/Indices ----------------------------- */

/* Get the active Wirecolor Set for the Scene provided */
WirecolorSet *SCENE_get_active_wirecolorset (Scene *scene)
{
	/* if no scene, we've got no hope of finding the Wirecolor Set */
	if (scene == NULL)
		return NULL;
	
	/* currently, there are two possibilities here:
	 *	-   0: no active Wirecolor set
	 *	- > 0: one of the user-defined Wirecolor Sets, but indices start from 0 (hence the -1)
	 */
	if (scene->active_wirecolorset > 0)
		return BLI_findlink(&scene->wirecolorsets, scene->active_wirecolorset-1);
	else
		return NULL;
}

/* Get the specified Wirecolor Set from the Scene provided */
WirecolorSet *SCENE_get_wirecolorset (Scene *scene, int idx)
{
	/* if no scene or zero index, we've got no hope of finding the Wirecolor Set */
	if (scene == NULL || idx == 0)
		return NULL;

	/* currently, there are two possibilities here:
	 *	-   0: no active Wirecolor set
	 *	- > 0: one of the user-defined Wirecolor Sets, but indices start from 0 (hence the -1)
	 */
	return BLI_findlink(&scene->wirecolorsets, idx-1);
}

/* Get the index of the Wirecolor Set provided, for the given Scene */
int SCENE_get_wirecolorset_index (Scene *scene, WirecolorSet *wcs)
{
	int index;
	
	/* if no WirecolorSet provided, have none */
	if (wcs == NULL)
		return 0;
	
	/* check if the WirecolorSet exists in scene list */
	if (scene) {
		/* get index and if valid, return 
		 *	- (absolute) Scene WirecolorSets are from (>= 1)
		 */
		index = BLI_findindex(&scene->wirecolorsets, wcs);
		if (index != -1)
			return (index + 1);
	}

	return 0;
}

/* Check if WirecolorSet can be used in the current context */
/* XXX - remove this? not needed */
short SCENE_wirecolorset_context_ok_poll (bContext *C, WirecolorSet *wcs)
{
	if ((wcs->flags & WCS_OVERRIDE_GROUP) == 0) {
		return 0;
	}

	return 1;
}

/* Menu of All Wirecolor Sets ----------------------------- */

/* Dynamically populate an enum of Wirecolor Sets */
EnumPropertyItem *SCENE_wirecolor_sets_enum_itemf (bContext *C, PointerRNA *UNUSED(ptr), PropertyRNA *UNUSED(prop), int *free)
{
	Scene *scene = CTX_data_scene(C);
	WirecolorSet *wcs;
	EnumPropertyItem *item= NULL, item_tmp= {0};
	int totitem= 0;
	int i= 0;

	if (C == NULL) {
		return DummyRNA_DEFAULT_items;
	}
	
	/* active Wirecolor Set
	 *	- only include entry if it exists
	 */
	if (scene->active_wirecolorset) {
		/* active Wirecolor Set */
		item_tmp.identifier= item_tmp.name= "Active Wirecolor Set";
		item_tmp.value= i++;
		RNA_enum_item_add(&item, &totitem, &item_tmp);
		
		/* separator */
		RNA_enum_item_add_separator(&item, &totitem);
	}
	else
		i++;
		
	/* user-defined Wirecolor Sets
	 *	- these are listed in the order in which they were defined for the active scene
	 */
	if (scene->wirecolorsets.first) {
		for (wcs= scene->wirecolorsets.first; wcs; wcs= wcs->next) {
			if (SCENE_wirecolorset_context_ok_poll(C, wcs)) {
				item_tmp.identifier= item_tmp.name= wcs->name;
				item_tmp.value= i++;
				RNA_enum_item_add(&item, &totitem, &item_tmp);
			}
		}
		
		/* separator */
		RNA_enum_item_add_separator(&item, &totitem);
	}
	
	RNA_enum_item_end(&item, &totitem);
	*free= 1;

	return item;
}

/* Create (and show) a menu containing all the Wirecolor Sets which can be used in the current context */
void SCENE_wirecolor_sets_menu_setup (bContext *C, const char title[], const char op_name[])
{
	Scene *scene= CTX_data_scene(C);
	WirecolorSet *wcs;
	uiPopupMenu *pup;
	uiLayout *layout;
	int i = 0;
	
	pup= uiPupMenuBegin(C, title, ICON_NONE);
	layout= uiPupMenuLayout(pup);
	
	/* active Wirecolor Set
	 *	- only include entry if it exists
	 */
	if (scene->active_wirecolorset) {
		uiItemIntO(layout, "Active Wirecolor Set", ICON_NONE, op_name, "type", i++);
		uiItemS(layout);
	}
	else
		i++;
	
	/* user-defined Wirecolor Sets
	 *	- these are listed in the order in which they were defined for the active scene
	 */
	if (scene->wirecolorsets.first) {
		for (wcs= scene->wirecolorsets.first; wcs; wcs=wcs->next, i++) {
			if (SCENE_wirecolorset_context_ok_poll(C, wcs))
				uiItemIntO(layout, wcs->name, ICON_NONE, op_name, "type", i);
		}
		uiItemS(layout);
	}
	
	uiPupMenuEnd(C, pup);
} 


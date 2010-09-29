/**
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
 * The Original Code is Copyright (C) 2009 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include <stdlib.h>
#include <string.h>

#include "MEM_guardedalloc.h"

#include "DNA_userdef_types.h"
#include "DNA_scene_types.h"

#include "BLI_fileops.h"
#include "BLI_path_util.h"
#include "BLI_storage.h"
#include "BLI_string.h"

#include "BKE_context.h"
#include "BKE_global.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_screen.h"

#include "RNA_access.h"

#include "UI_interface.h"

#include "buttons_intern.h"	// own include

/********************** toolbox operator *********************/

static int toolbox_invoke(bContext *C, wmOperator *op, wmEvent *event)
{
	bScreen *sc= CTX_wm_screen(C);
	SpaceButs *sbuts= CTX_wm_space_buts(C);
	PointerRNA ptr;
	uiPopupMenu *pup;
	uiLayout *layout;

	RNA_pointer_create(&sc->id, &RNA_SpaceProperties, sbuts, &ptr);

	pup= uiPupMenuBegin(C, "Align", 0);
	layout= uiPupMenuLayout(pup);
	uiItemsEnumR(layout, &ptr, "align");
	uiPupMenuEnd(C, pup);

	return OPERATOR_CANCELLED;
}

void BUTTONS_OT_toolbox(wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Toolbox";
	ot->description="Display button panel toolbox";
	ot->idname= "BUTTONS_OT_toolbox";
	
	/* api callbacks */
	ot->invoke= toolbox_invoke;
	ot->poll= ED_operator_buttons_active;
}

/********************** resflip operator *********************/

static int resflip_invoke(bContext *C, wmOperator *op, wmEvent *event)
{
	Scene *sce= CTX_data_scene(C);
	int tmpres;
	
	tmpres= sce->r.xsch;
	sce->r.xsch= sce->r.ysch;
	sce->r.ysch= tmpres;
	
	WM_event_add_notifier(C, NC_SCENE|ND_SPACE_PROPERTIES, CTX_data_active_object(C));
	
	return OPERATOR_FINISHED;
}

void BUTTONS_OT_resflip(wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Resolution Flip";
	ot->description="Swap horiz and vert resolutions - switch portrait/landscape.";
	ot->idname= "BUTTONS_OT_resflip";
	
	/* api callbacks */
	ot->invoke= resflip_invoke;
	ot->poll= ED_operator_buttons_active;
}

/********************** filebrowse operator *********************/

typedef struct FileBrowseOp {
	PointerRNA ptr;
	PropertyRNA *prop;
} FileBrowseOp;

static int file_browse_exec(bContext *C, wmOperator *op)
{
	FileBrowseOp *fbo= op->customdata;
	ID *id;
	char *base, *str, path[FILE_MAX];
	
	if (RNA_property_is_set(op->ptr, "filepath")==0 || fbo==NULL)
		return OPERATOR_CANCELLED;
	
	str= RNA_string_get_alloc(op->ptr, "filepath", 0, 0);

	/* add slash for directories, important for some properties */
	if(RNA_property_subtype(fbo->prop) == PROP_DIRPATH) {
		id = fbo->ptr.id.data;
		base = (id && id->lib)? id->lib->filepath: G.sce;

		BLI_strncpy(path, str, FILE_MAX);
		BLI_path_abs(path, base);

		if(BLI_is_dir(path)) {
			str = MEM_reallocN(str, strlen(str)+2);
			BLI_add_slash(str);
		}
	}

	RNA_property_string_set(&fbo->ptr, fbo->prop, str);
	RNA_property_update(C, &fbo->ptr, fbo->prop);
	MEM_freeN(str);

	MEM_freeN(op->customdata);
	return OPERATOR_FINISHED;
}

static int file_browse_cancel(bContext *C, wmOperator *op)
{
	MEM_freeN(op->customdata);
	op->customdata= NULL;

	return OPERATOR_CANCELLED;
}

static int file_browse_invoke(bContext *C, wmOperator *op, wmEvent *event)
{
	PointerRNA ptr;
	PropertyRNA *prop;
	FileBrowseOp *fbo;
	char *str;

	uiFileBrowseContextProperty(C, &ptr, &prop);

	if(!prop)
		return OPERATOR_CANCELLED;

	str= RNA_property_string_get_alloc(&ptr, prop, 0, 0);

	/* useful yet irritating feature, Shift+Click to open the file
	 * Alt+Click to browse a folder in the OS's browser */
	if(event->shift || event->alt) {
		PointerRNA props_ptr;

		if(event->alt) {
			char *lslash= BLI_last_slash(str);
			if(lslash)
				*lslash= '\0';
		}


		WM_operator_properties_create(&props_ptr, "WM_OT_path_open");
		RNA_string_set(&props_ptr, "filepath", str);
		WM_operator_name_call(C, "WM_OT_path_open", WM_OP_EXEC_DEFAULT, &props_ptr);
		WM_operator_properties_free(&props_ptr);

		MEM_freeN(str);
		return OPERATOR_CANCELLED;
	}
	else {
		fbo= MEM_callocN(sizeof(FileBrowseOp), "FileBrowseOp");
		fbo->ptr= ptr;
		fbo->prop= prop;
		op->customdata= fbo;

		RNA_string_set(op->ptr, "filepath", str);
		MEM_freeN(str);

		if(RNA_struct_find_property(op->ptr, "relative_path"))
			if(!RNA_property_is_set(op->ptr, "relative_path"))
				RNA_boolean_set(op->ptr, "relative_path", U.flag & USER_RELPATHS);

		WM_event_add_fileselect(C, op);

		return OPERATOR_RUNNING_MODAL;
	}
}

void BUTTONS_OT_file_browse(wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Accept";
	ot->description="Open a file browser, Hold Shift to open the file, Alt to browse containing directory";
	ot->idname= "BUTTONS_OT_file_browse";
	
	/* api callbacks */
	ot->invoke= file_browse_invoke;
	ot->exec= file_browse_exec;
	ot->cancel= file_browse_cancel;

	/* properties */
	WM_operator_properties_filesel(ot, 0, FILE_SPECIAL, FILE_OPENFILE, WM_FILESEL_FILEPATH|WM_FILESEL_RELPATH);
}


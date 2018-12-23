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
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/blenloader/intern/versioning_userdef.c
 *  \ingroup blenloader
 *
 * Version patch user preferences.
 */

#include <string.h>

#include "BLI_math.h"
#include "BLI_utildefines.h"

#include "DNA_userdef_types.h"
#include "DNA_curve_types.h"
#include "DNA_windowmanager_types.h"

#include "BKE_addon.h"
#include "BKE_colorband.h"
#include "BKE_idprop.h"
#include "BKE_main.h"
#include "BKE_keyconfig.h"

#include "BLO_readfile.h"  /* Own include. */

#include "wm_event_types.h"

/* Disallow access to global userdef. */
#define U (_error_)


static void do_versions_theme(UserDef *userdef, bTheme *btheme)
{

#define USER_VERSION_ATLEAST(ver, subver) MAIN_VERSION_ATLEAST(userdef, ver, subver)
	if (!USER_VERSION_ATLEAST(280, 20)) {
		memcpy(btheme, &U_theme_default, sizeof(*btheme));
	}

	if (!USER_VERSION_ATLEAST(280, 25)) {
		copy_v4_v4_char(btheme->tact.anim_preview_range, btheme->tact.anim_active);
		copy_v4_v4_char(btheme->tnla.anim_preview_range, btheme->tnla.anim_active);
		copy_v4_v4_char(btheme->tipo.anim_preview_range, btheme->tact.anim_active);
	}

	if (!USER_VERSION_ATLEAST(280, 26)) {
		copy_v4_v4_char(btheme->tui.icon_collection, U_theme_default.tui.icon_collection);
		copy_v4_v4_char(btheme->tui.icon_object, U_theme_default.tui.icon_object);
		copy_v4_v4_char(btheme->tui.icon_object_data, U_theme_default.tui.icon_object_data);
		copy_v4_v4_char(btheme->tui.icon_modifier, U_theme_default.tui.icon_modifier);
		copy_v4_v4_char(btheme->tui.icon_shading, U_theme_default.tui.icon_shading);
	}

	if (!USER_VERSION_ATLEAST(280, 27)) {
		copy_v4_v4_char(btheme->tact.shade2, U_theme_default.tact.shade2);
		copy_v4_v4_char(btheme->tact.hilite, U_theme_default.tact.hilite);
		copy_v4_v4_char(btheme->tact.group, U_theme_default.tact.group);
		copy_v4_v4_char(btheme->tact.group_active, U_theme_default.tact.group_active);
		copy_v4_v4_char(btheme->tact.strip_select, U_theme_default.tact.strip_select);
		copy_v4_v4_char(btheme->tact.ds_channel, U_theme_default.tact.ds_channel);
		copy_v4_v4_char(btheme->tact.ds_subchannel, U_theme_default.tact.ds_subchannel);
		copy_v4_v4_char(btheme->tact.keytype_movehold, U_theme_default.tact.keytype_movehold);
		copy_v4_v4_char(btheme->tact.keytype_movehold_select, U_theme_default.tact.keytype_movehold_select);
	}

	if (!USER_VERSION_ATLEAST(280, 28)) {
		copy_v4_v4_char(btheme->tact.ds_ipoline, U_theme_default.tact.ds_ipoline);
	}

	if (!USER_VERSION_ATLEAST(280, 29)) {
		copy_v4_v4_char(btheme->tbuts.navigation_bar, U_theme_default.ttopbar.header);
	}
	if (!USER_VERSION_ATLEAST(280, 31)) {
		copy_v4_v4_char(btheme->tclip.list_text, U_theme_default.tclip.list_text);
	}

	if (!USER_VERSION_ATLEAST(280, 33)) {
		copy_v4_v4_char(btheme->tuserpref.navigation_bar, U_theme_default.tuserpref.navigation_bar);
	}

	if (!USER_VERSION_ATLEAST(280, 36)) {
		copy_v4_v4_char(btheme->tui.wcol_state.inner_changed, U_theme_default.tui.wcol_state.inner_changed);
		copy_v4_v4_char(btheme->tui.wcol_state.inner_changed_sel, U_theme_default.tui.wcol_state.inner_changed_sel);
	}

	if (!USER_VERSION_ATLEAST(280, 39)) {
		copy_v4_v4_char(btheme->tclip.metadatabg, U_theme_default.tima.metadatabg);
		copy_v4_v4_char(btheme->tclip.metadatatext, U_theme_default.tima.metadatatext);
	}
#undef USER_VERSION_ATLEAST
}

/* UserDef.flag */
#define USER_LMOUSESELECT (1 << 14)  /* deprecated */

static void do_version_select_mouse(UserDef *userdef, wmKeyMapItem *kmi)
{
	/* Remove select/action mouse from user defined keymaps. */
	enum {
		ACTIONMOUSE = 0x0005,
		SELECTMOUSE = 0x0006,
		EVT_TWEAK_A = 0x5005,
		EVT_TWEAK_S = 0x5006,
	};
	const bool left = (userdef->flag & USER_LMOUSESELECT) != 0;

	switch (kmi->type) {
		case SELECTMOUSE: kmi->type = (left) ? LEFTMOUSE : RIGHTMOUSE; break;
		case ACTIONMOUSE: kmi->type = (left) ? RIGHTMOUSE : LEFTMOUSE; break;
		case EVT_TWEAK_S: kmi->type = (left) ? EVT_TWEAK_L : EVT_TWEAK_R; break;
		case EVT_TWEAK_A: kmi->type = (left) ? EVT_TWEAK_R : EVT_TWEAK_L; break;
		default: break;
	}
}

/* patching UserDef struct and Themes */
void BLO_version_defaults_userpref_blend(Main *bmain, UserDef *userdef)
{

#define USER_VERSION_ATLEAST(ver, subver) MAIN_VERSION_ATLEAST(bmain, ver, subver)

	/* the UserDef struct is not corrected with do_versions() .... ugh! */
	if (userdef->wheellinescroll == 0) userdef->wheellinescroll = 3;
	if (userdef->menuthreshold1 == 0) {
		userdef->menuthreshold1 = 5;
		userdef->menuthreshold2 = 2;
	}
	if (userdef->tb_leftmouse == 0) {
		userdef->tb_leftmouse = 5;
		userdef->tb_rightmouse = 5;
	}
	if (userdef->mixbufsize == 0) userdef->mixbufsize = 2048;
	if (userdef->autokey_mode == 0) {
		/* 'add/replace' but not on */
		userdef->autokey_mode = 2;
	}
	if (userdef->savetime <= 0) {
		userdef->savetime = 1;
// XXX		error(STRINGIFY(BLENDER_STARTUP_FILE)" is buggy, please consider removing it.\n");
	}
	if (userdef->gizmo_size == 0) {
		userdef->gizmo_size = 75;
		userdef->gizmo_flag |= USER_GIZMO_DRAW;
	}
	if (userdef->pad_rot_angle == 0.0f)
		userdef->pad_rot_angle = 15.0f;

	/* graph editor - unselected F-Curve visibility */
	if (userdef->fcu_inactive_alpha == 0) {
		userdef->fcu_inactive_alpha = 0.25f;
	}

	if (!USER_VERSION_ATLEAST(192, 0)) {
		strcpy(userdef->sounddir, "/");
	}

	/* patch to set Dupli Armature */
	if (!USER_VERSION_ATLEAST(220, 0)) {
		userdef->dupflag |= USER_DUP_ARM;
	}

	/* added seam, normal color, undo */
	if (!USER_VERSION_ATLEAST(235, 0)) {
		userdef->uiflag |= USER_GLOBALUNDO;
		if (userdef->undosteps == 0) userdef->undosteps = 32;
	}
	if (!USER_VERSION_ATLEAST(236, 0)) {
		/* illegal combo... */
		if (userdef->flag & USER_LMOUSESELECT)
			userdef->flag &= ~USER_TWOBUTTONMOUSE;
	}
	if (!USER_VERSION_ATLEAST(240, 0)) {
		userdef->uiflag |= USER_PLAINMENUS;
		if (userdef->obcenter_dia == 0) userdef->obcenter_dia = 6;
	}
	if (!USER_VERSION_ATLEAST(242, 0)) {
		/* set defaults for 3D View rotating axis indicator */
		/* since size can't be set to 0, this indicates it's not saved in startup.blend */
		if (userdef->rvisize == 0) {
			userdef->rvisize = 15;
			userdef->rvibright = 8;
			userdef->uiflag |= USER_SHOW_GIZMO_AXIS;
		}

	}
	if (!USER_VERSION_ATLEAST(244, 0)) {
		/* set default number of recently-used files (if not set) */
		if (userdef->recent_files == 0) userdef->recent_files = 10;
	}
	if (!USER_VERSION_ATLEAST(245, 3)) {
		if (userdef->coba_weight.tot == 0)
			BKE_colorband_init(&userdef->coba_weight, true);
	}
	if (!USER_VERSION_ATLEAST(245, 3)) {
		userdef->flag |= USER_ADD_VIEWALIGNED | USER_ADD_EDITMODE;
	}
	if (!USER_VERSION_ATLEAST(250, 0)) {
		/* adjust grease-pencil distances */
		userdef->gp_manhattendist = 1;
		userdef->gp_euclideandist = 2;

		/* adjust default interpolation for new IPO-curves */
		userdef->ipo_new = BEZT_IPO_BEZ;
	}

	if (!USER_VERSION_ATLEAST(250, 3)) {
		/* new audio system */
		if (userdef->audiochannels == 0)
			userdef->audiochannels = 2;
		if (userdef->audioformat == 0)
			userdef->audioformat = 0x24;
		if (userdef->audiorate == 0)
			userdef->audiorate = 48000;
	}

	if (!USER_VERSION_ATLEAST(250, 8)) {
		wmKeyMap *km;

		for (km = userdef->user_keymaps.first; km; km = km->next) {
			if (STREQ(km->idname, "Armature_Sketch"))
				strcpy(km->idname, "Armature Sketch");
			else if (STREQ(km->idname, "View3D"))
				strcpy(km->idname, "3D View");
			else if (STREQ(km->idname, "View3D Generic"))
				strcpy(km->idname, "3D View Generic");
			else if (STREQ(km->idname, "EditMesh"))
				strcpy(km->idname, "Mesh");
			else if (STREQ(km->idname, "UVEdit"))
				strcpy(km->idname, "UV Editor");
			else if (STREQ(km->idname, "Animation_Channels"))
				strcpy(km->idname, "Animation Channels");
			else if (STREQ(km->idname, "GraphEdit Keys"))
				strcpy(km->idname, "Graph Editor");
			else if (STREQ(km->idname, "GraphEdit Generic"))
				strcpy(km->idname, "Graph Editor Generic");
			else if (STREQ(km->idname, "Action_Keys"))
				strcpy(km->idname, "Dopesheet");
			else if (STREQ(km->idname, "NLA Data"))
				strcpy(km->idname, "NLA Editor");
			else if (STREQ(km->idname, "Node Generic"))
				strcpy(km->idname, "Node Editor");
			else if (STREQ(km->idname, "Logic Generic"))
				strcpy(km->idname, "Logic Editor");
			else if (STREQ(km->idname, "File"))
				strcpy(km->idname, "File Browser");
			else if (STREQ(km->idname, "FileMain"))
				strcpy(km->idname, "File Browser Main");
			else if (STREQ(km->idname, "FileButtons"))
				strcpy(km->idname, "File Browser Buttons");
			else if (STREQ(km->idname, "Buttons Generic"))
				strcpy(km->idname, "Property Editor");
		}
	}

	if (!USER_VERSION_ATLEAST(252, 3)) {
		if (userdef->flag & USER_LMOUSESELECT)
			userdef->flag &= ~USER_TWOBUTTONMOUSE;
	}
	if (!USER_VERSION_ATLEAST(252, 4)) {
		/* default new handle type is auto handles */
		userdef->keyhandles_new = HD_AUTO;
	}

	if (!USER_VERSION_ATLEAST(257, 0)) {
		/* clear "AUTOKEY_FLAG_ONLYKEYINGSET" flag from userprefs,
		 * so that it doesn't linger around from old configs like a ghost */
		userdef->autokey_flag &= ~AUTOKEY_FLAG_ONLYKEYINGSET;
	}

	if (!USER_VERSION_ATLEAST(260, 3)) {
		/* if new keyframes handle default is stuff "auto", make it "auto-clamped" instead
		 * was changed in 260 as part of GSoC11, but version patch was wrong
		 */
		if (userdef->keyhandles_new == HD_AUTO)
			userdef->keyhandles_new = HD_AUTO_ANIM;

		/* enable (Cycles) addon by default */
		BKE_addon_ensure(&userdef->addons, "cycles");
	}

	if (!USER_VERSION_ATLEAST(261, 4)) {
		userdef->use_16bit_textures = true;
	}

	if (!USER_VERSION_ATLEAST(267, 0)) {

		/* GL Texture Garbage Collection */
		if (userdef->textimeout == 0) {
			userdef->texcollectrate = 60;
			userdef->textimeout = 120;
		}
		if (userdef->memcachelimit <= 0) {
			userdef->memcachelimit = 32;
		}
		if (userdef->dbl_click_time == 0) {
			userdef->dbl_click_time = 350;
		}
		if (userdef->v2d_min_gridsize == 0) {
			userdef->v2d_min_gridsize = 35;
		}
		if (userdef->dragthreshold == 0)
			userdef->dragthreshold = 5;
		if (userdef->widget_unit == 0)
			userdef->widget_unit = 20;
		if (userdef->anisotropic_filter <= 0)
			userdef->anisotropic_filter = 1;

		if (userdef->ndof_sensitivity == 0.0f) {
			userdef->ndof_sensitivity = 1.0f;
			userdef->ndof_flag = (NDOF_LOCK_HORIZON | NDOF_SHOULD_PAN | NDOF_SHOULD_ZOOM | NDOF_SHOULD_ROTATE);
		}

		if (userdef->ndof_orbit_sensitivity == 0.0f) {
			userdef->ndof_orbit_sensitivity = userdef->ndof_sensitivity;

			if (!(userdef->flag & USER_TRACKBALL))
				userdef->ndof_flag |= NDOF_TURNTABLE;
		}
		if (userdef->tweak_threshold == 0)
			userdef->tweak_threshold = 10;
	}

	/* NOTE!! from now on use userdef->versionfile and userdef->subversionfile */
#undef USER_VERSION_ATLEAST
#define USER_VERSION_ATLEAST(ver, subver) MAIN_VERSION_ATLEAST(userdef, ver, subver)

	if (!USER_VERSION_ATLEAST(271, 5)) {
		userdef->pie_menu_radius = 100;
		userdef->pie_menu_threshold = 12;
		userdef->pie_animation_timeout = 6;
	}

	if (!USER_VERSION_ATLEAST(275, 2)) {
		userdef->ndof_deadzone = 0.1;
	}

	if (!USER_VERSION_ATLEAST(275, 4)) {
		userdef->node_margin = 80;
	}

	if (!USER_VERSION_ATLEAST(278, 6)) {
		/* Clear preference flags for re-use. */
		userdef->flag &= ~(
		    USER_FLAG_NUMINPUT_ADVANCED | USER_FLAG_DEPRECATED_2 | USER_FLAG_DEPRECATED_3 |
		    USER_FLAG_DEPRECATED_6 | USER_FLAG_DEPRECATED_7 |
		    USER_FLAG_DEPRECATED_9 | USER_DEVELOPER_UI);
		userdef->uiflag &= ~(
		    USER_HEADER_BOTTOM);
		userdef->transopts &= ~(
		    USER_TR_DEPRECATED_2 | USER_TR_DEPRECATED_3 | USER_TR_DEPRECATED_4 |
		    USER_TR_DEPRECATED_6 | USER_TR_DEPRECATED_7);

		userdef->uiflag |= USER_LOCK_CURSOR_ADJUST;
	}


	if (!USER_VERSION_ATLEAST(280, 20)) {
		userdef->gpu_viewport_quality = 0.6f;

		/* Reset theme, old themes will not be compatible with minor version updates from now on. */
		for (bTheme *btheme = userdef->themes.first; btheme; btheme = btheme->next) {
			memcpy(btheme, &U_theme_default, sizeof(*btheme));
		}

		/* Annotations - new layer color
		 * Replace anything that used to be set if it looks like was left
		 * on the old default (i.e. black), which most users used
		 */
		if ((userdef->gpencil_new_layer_col[3] < 0.1f) || (userdef->gpencil_new_layer_col[0] < 0.1f)) {
			/* - New color matches the annotation pencil icon
			 * - Non-full alpha looks better!
			 */
			ARRAY_SET_ITEMS(userdef->gpencil_new_layer_col, 0.38f, 0.61f, 0.78f, 0.9f);
		}
	}

	if (!USER_VERSION_ATLEAST(280, 31)) {
		/* Remove select/action mouse from user defined keymaps. */
		for (wmKeyMap *keymap = userdef->user_keymaps.first; keymap; keymap = keymap->next) {
			for (wmKeyMapDiffItem *kmdi = keymap->diff_items.first; kmdi; kmdi = kmdi->next) {
				if (kmdi->remove_item) {
					do_version_select_mouse(userdef, kmdi->remove_item);
				}
				if (kmdi->add_item) {
					do_version_select_mouse(userdef, kmdi->add_item);
				}
			}

			for (wmKeyMapItem *kmi = keymap->items.first; kmi; kmi = kmi->next) {
				do_version_select_mouse(userdef, kmi);
			}
		}
	}

	if (!USER_VERSION_ATLEAST(280, 33)) {
		/* Enable GLTF addon by default. */
		BKE_addon_ensure(&userdef->addons, "io_scene_gltf2");
	}

	if (!USER_VERSION_ATLEAST(280, 35)) {
		/* Preserve RMB select setting after moving to Python and changing default value. */
		if (USER_VERSION_ATLEAST(280, 32) || !(userdef->flag & USER_LMOUSESELECT)) {
			BKE_keyconfig_pref_set_select_mouse(userdef, 1, false);
		}

		userdef->flag &= ~USER_LMOUSESELECT;
	}

	if (!USER_VERSION_ATLEAST(280, 38)) {

		/* (keep this block even if it becomes empty). */
		copy_v4_fl4(userdef->light_param[0].vec, -0.580952, 0.228571, 0.781185, 0.0);
		copy_v4_fl4(userdef->light_param[0].col, 0.900000, 0.900000, 0.900000, 1.000000);
		copy_v4_fl4(userdef->light_param[0].spec, 0.318547, 0.318547, 0.318547, 1.000000);
		userdef->light_param[0].flag = 1;
		userdef->light_param[0].smooth = 0.1;

		copy_v4_fl4(userdef->light_param[1].vec, 0.788218, 0.593482, -0.162765, 0.0);
		copy_v4_fl4(userdef->light_param[1].col, 0.267115, 0.269928, 0.358840, 1.000000);
		copy_v4_fl4(userdef->light_param[1].spec, 0.090838, 0.090838, 0.090838, 1.000000);
		userdef->light_param[1].flag = 1;
		userdef->light_param[1].smooth = 0.25;

		copy_v4_fl4(userdef->light_param[2].vec, 0.696472, -0.696472, -0.172785, 0.0);
		copy_v4_fl4(userdef->light_param[2].col, 0.293216, 0.304662, 0.401968, 1.000000);
		copy_v4_fl4(userdef->light_param[2].spec, 0.069399, 0.020331, 0.020331, 1.000000);
		userdef->light_param[2].flag = 1;
		userdef->light_param[2].smooth = 0.4;

		copy_v4_fl4(userdef->light_param[3].vec, 0.021053, -0.989474, 0.143173, 0.0);
		copy_v4_fl4(userdef->light_param[3].col, 0.0, 0.0, 0.0, 1.0);
		copy_v4_fl4(userdef->light_param[3].spec, 0.072234, 0.082253, 0.162642, 1.000000);
		userdef->light_param[3].flag = 1;
		userdef->light_param[3].smooth = 0.7;

		copy_v4_fl4(userdef->light_ambient, 0.025000, 0.025000, 0.025000, 1.000000);

		userdef->flag &= ~(
		        USER_FLAG_DEPRECATED_4);

		userdef->uiflag &= ~(
		        USER_UIFLAG_DEPRECATED_8 |
		        USER_UIFLAG_DEPRECATED_12 |
		        USER_UIFLAG_DEPRECATED_22);
	}

	/**
	 * Include next version bump.
	 */
	{
		/* (keep this block even if it becomes empty). */
	}

	if (userdef->pixelsize == 0.0f)
		userdef->pixelsize = 1.0f;

	if (userdef->image_draw_method == 0)
		userdef->image_draw_method = IMAGE_DRAW_METHOD_2DTEXTURE;

	// we default to the first audio device
	userdef->audiodevice = 0;

	for (bTheme *btheme = userdef->themes.first; btheme; btheme = btheme->next) {
		do_versions_theme(userdef, btheme);
	}
#undef USER_VERSION_ATLEAST

}

#undef USER_LMOUSESELECT

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
 * The Original Code is Copyright (C) 2013, Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file ED_buttons.h
 *  \ingroup editors
 */

#ifndef __ED_BUTTONS_H__
#define __ED_BUTTONS_H__

#include "BLI_utildefines.h"

/* Used to check whether a given texture context is valid in current context. */
bool ED_texture_context_check_world(const struct bContext *C);
bool ED_texture_context_check_material(const struct bContext *C);
bool ED_texture_context_check_lamp(const struct bContext *C);
bool ED_texture_context_check_particles(const struct bContext *C);
bool ED_texture_context_check_linestyle(const struct bContext *C);
bool ED_texture_context_check_others(const struct bContext *C);

void ED_buttons_id_unref(struct SpaceButs *sbuts, const struct ID *id);

#endif /*  __ED_BUTTONS_H__ */

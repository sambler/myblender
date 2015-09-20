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
 * The Original Code is Copyright (C) 2015, Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file ED_outliner.h
 *  \ingroup editors
 */

#ifndef __ED_OUTLINER_H__
#define __ED_OUTLINER_H__

struct ID;
struct SpaceOops;

/* Used to check whether a given texture context is valid in current context. */
void ED_outliner_id_unref(struct SpaceOops *so, const struct ID *id);

#endif /*  __ED_OUTLINER_H__ */

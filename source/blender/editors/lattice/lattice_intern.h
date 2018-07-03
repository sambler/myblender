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
 * The Original Code is Copyright (C) 2008 Blender Foundation.
 * All rights reserved.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/lattice/lattice_intern.h
 *  \ingroup edlattice
 */


#ifndef __LATTICE_INTERN_H__
#define __LATTICE_INTERN_H__

/* editlattice_select.c */
void LATTICE_OT_select_all(struct wmOperatorType *ot);
void LATTICE_OT_select_more(struct wmOperatorType *ot);
void LATTICE_OT_select_less(struct wmOperatorType *ot);
void LATTICE_OT_select_ungrouped(struct wmOperatorType *ot);
void LATTICE_OT_select_random(struct wmOperatorType *ot);
void LATTICE_OT_select_mirror(struct wmOperatorType *ot);

/* editlattice_tools.c */
void LATTICE_OT_make_regular(struct wmOperatorType *ot);
void LATTICE_OT_flip(struct wmOperatorType *ot);

#endif  /* __LATTICE_INTERN_H__ */

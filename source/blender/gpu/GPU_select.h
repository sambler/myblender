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
 * The Original Code is Copyright (C) 2014 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Antony Riakiotakis.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file GPU_select.h
 *  \ingroup gpu
 */

#ifndef __GPU_SELECT_H__
#define __GPU_SELECT_H__

#include "DNA_vec_types.h"  /* rcft */
#include "BLI_sys_types.h"

/* flags for mode of operation */
enum {
	GPU_SELECT_ALL                      = 1,
	GPU_SELECT_NEAREST_FIRST_PASS       = 2,
	GPU_SELECT_NEAREST_SECOND_PASS      = 3,
};

void GPU_select_begin(unsigned int *buffer, unsigned int bufsize, rctf *input, char mode, int oldhits);
bool GPU_select_load_id(unsigned int id);
unsigned int GPU_select_end(void);
bool GPU_select_query_check_active(void);

#endif

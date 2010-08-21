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
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 * Unix specific implementations
 */

#ifndef KXH_UNIX_SERVICES_H
#define KXH_UNIX_SERVICES_H

#include "KXH_engine_data_wraps.h"

#ifdef __cplusplus
extern "C" {
#endif

	/** Create appropriate devices for the gameengine. */
	void
	KXH_create_devices(
		ketsji_engine_data* k
		);

	/** Add banners to the canvas. */
	void
	KXH_add_banners(
		ketsji_engine_data* k
		);


#ifdef __cplusplus
}
#endif

#endif


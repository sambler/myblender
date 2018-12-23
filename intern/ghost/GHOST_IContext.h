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
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file ghost/GHOST_IContext.h
 *  \ingroup GHOST
 * Declaration of GHOST_IContext interface class.
 */

#ifndef __GHOST_IContext_H__
#define __GHOST_IContext_H__

#include "STR_String.h"
#include "GHOST_Types.h"


/**
 * Interface for GHOST context.
 *
 * You can create a offscreen context (windowless) with the system's
 * GHOST_ISystem::createOffscreenContext method.
 * \see GHOST_ISystem#createOffscreenContext
 *
 * \author  Clément Foucault
 * \date    Feb 9, 2018
 */
class GHOST_IContext
{
public:
	/**
	 * Destructor.
	 */
	virtual ~GHOST_IContext()
	{
	}

	/**
	 * Activates the drawing context.
	 * \return  A boolean success indicator.
	 */
	virtual GHOST_TSuccess activateDrawingContext() = 0;

	/**
	 * Release the drawing context of the calling thread.
	 * \return  A boolean success indicator.
	 */
	virtual GHOST_TSuccess releaseDrawingContext() = 0;

#ifdef WITH_CXX_GUARDEDALLOC
	MEM_CXX_CLASS_ALLOC_FUNCS("GHOST:GHOST_IContext")
#endif
};

#endif // __GHOST_IContext_H__

# $Id$
# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# The Original Code is Copyright (C) 2011, Blender Foundation
# All rights reserved.
#
# The Original Code is: all of this file.
#
# Contributor(s): Shane Ambler.
#
# ***** END GPL LICENSE BLOCK *****

set(INC 
	.
	../blenkernel
	../blenlib
	../blenloader
	../imbuf
	../makesdna
	../makesrna
	../../../intern/guardedalloc
)

set(INC_SYS

)

set(SRC
	intern/add.c
	intern/adjustment.c
	intern/alpha.c
	intern/cross.c
	intern/drop.c
	intern/gamma.c
	intern/glow.c
	intern/mul.c
	intern/multicam.c
	intern/plugin.c
	intern/seqeffects.c
	intern/solid_colour.c
	intern/speed.c
	intern/sub.c
	intern/transform.c
	intern/utils.c
	intern/wipe.c
	
	intern/seq_intern.h
	SEQ_effects.h
)

if(WITH_OPENMP)
	add_definitions(-DPARALLEL=1)
endif()

blender_add_lib(bf_sequencer_effects "${SRC}" "${INC}" "${INC_SYS}")

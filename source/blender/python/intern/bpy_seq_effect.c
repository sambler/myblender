
/*
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
 * Contributor(s): Shane Ambler
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/python/intern/bpy_seq_effect.c
 *  \ingroup pythonintern
 */


/* Note, this module is not to be used directly by the user.
 * Internally its exposed as '_bpy.seqfx', which provides functions for 'bpy.seqfx', a python package.
 * */

#include <Python.h>

#include "RNA_types.h"

#include "BPY_extern.h"
#include "bpy_seq_effect.h"
#include "bpy_rna.h" /* for setting arg props only - pyrna_py_to_prop() */
#include "bpy_util.h"
#include "../generic/bpy_internal_import.h"

#include "BLI_utildefines.h"

#include "RNA_access.h"
#include "RNA_enum_types.h"

#include "WM_api.h"
#include "WM_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_ghash.h"

#include "BKE_report.h"
#include "BKE_context.h"


static PyObject *pyseqfx_init(PyObject *UNUSED(self), PyObject *args)
{
	
}

static PyObject *pyseqfx_num_inputs(PyObject *UNUSED(self), PyObject *args)
{
	
}

static PyObject *pyseqfx_load(PyObject *UNUSED(self), PyObject *args)
{
	
}

static PyObject *pyseqfx_copy(PyObject *UNUSED(self))
{
	
}

static PyObject *pyseqfx_free(PyObject *UNUSED(self), PyObject *value)
{
	
}

static PyObject *pyseqfx_early_out(PyObject *UNUSED(self), PyObject *value)
{

}

static PyObject *pyseqfx_store_icu_range(PyObject *UNUSED(self), PyObject *value)
{

}

static PyObject *pyseqfx_get_default_fac(PyObject *UNUSED(self), PyObject *value)
{

}

static PyObject *pyseqfx_execute(PyObject *UNUSED(self), PyObject *value)
{

}

static struct PyMethodDef bpy_seqfx_methods[]= {
	{"init", (PyCFunction) pyseqfx_init, METH_VARARGS, NULL},
	{"num_inputs", (PyCFunction) pyseqfx_num_inputs, METH_VARARGS, NULL},
	{"load", (PyCFunction) pyseqfx_load, METH_VARARGS, NULL},
	{"copy", (PyCFunction) pyseqfx_copy, METH_VARARGS, NULL},
	{"free", (PyCFunction) pyseqfx_free, METH_VARARGS, NULL},
	{"early_out", (PyCFunction) pyseqfx_early_out, METH_VARARGS, NULL},
	{"store_icu_range", (PyCFunction) pyseqfx_store_icu_range, METH_VARARGS, NULL},
	{"get_default_fac", (PyCFunction) pyseqfx_get_default_fac, METH_VARARGS, NULL},
	{"execute", (PyCFunction) pyseqfx_execute, METH_VARARGS, NULL},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef bpy_seqfx_module= {
	PyModuleDef_HEAD_INIT,
	"_bpy.seqfx",
	NULL,
	-1,/* multiple "initialization" just copies the module dict. */
	bpy_seqfx_methods,
	NULL, NULL, NULL, NULL
};

PyObject *BPY_seqfx_module(void)
{
	PyObject *submodule;

	submodule= PyModule_Create(&bpy_seqfx_module);

	return submodule;
}

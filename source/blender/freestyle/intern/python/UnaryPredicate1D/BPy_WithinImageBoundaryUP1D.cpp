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
 */

/** \file source/blender/freestyle/intern/python/UnaryPredicate1D/BPy_WithinImageBoundaryUP1D.cpp
 *  \ingroup freestyle
 */

#include "BPy_WithinImageBoundaryUP1D.h"

#ifdef __cplusplus
extern "C" {
#endif

///////////////////////////////////////////////////////////////////////////////////////////

//------------------------INSTANCE METHODS ----------------------------------

static char WithinImageBoundaryUP1D___doc__[] =
"Class hierarchy: :class:`freestyle.types.UnaryPredicate1D` > :class:`WithinImageBoundaryUP1D`\n"
"\n"
".. method:: __init__(xmin, ymin, xmax, ymax)\n"
"\n"
"   Builds an WithinImageBoundaryUP1D object.\n"
"\n"
"   :arg xmin: X lower bound of the image boundary.\n"
"   :type xmin: float\n"
"   :arg ymin: Y lower bound of the image boundary.\n"
"   :type ymin: float\n"
"   :arg xmax: X upper bound of the image boundary.\n"
"   :type xmax: float\n"
"   :arg ymax: Y upper bound of the image boundary.\n"
"   :type ymax: float\n"
"\n"
".. method:: __call__(inter)\n"
"\n"
"   Returns true if the Interface1D intersects with image boundary.\n";

static int WithinImageBoundaryUP1D___init__(BPy_WithinImageBoundaryUP1D *self, PyObject *args, PyObject *kwds)
{
	static const char *kwlist[] = {"xmin", "ymin", "xmax", "ymax", NULL};
	double xmin, ymin, xmax, ymax;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddd", (char **)kwlist, &xmin, &ymin, &xmax, &ymax))
		return -1;
	self->py_up1D.up1D = new Predicates1D::WithinImageBoundaryUP1D(xmin, ymin, xmax, ymax);
	return 0;
}

/*-----------------------BPy_TrueUP1D type definition ------------------------------*/

PyTypeObject WithinImageBoundaryUP1D_Type = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"WithinImageBoundaryUP1D",      /* tp_name */
	sizeof(BPy_WithinImageBoundaryUP1D), /* tp_basicsize */
	0,                              /* tp_itemsize */
	0,                              /* tp_dealloc */
	0,                              /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_reserved */
	0,                              /* tp_repr */
	0,                              /* tp_as_number */
	0,                              /* tp_as_sequence */
	0,                              /* tp_as_mapping */
	0,                              /* tp_hash  */
	0,                              /* tp_call */
	0,                              /* tp_str */
	0,                              /* tp_getattro */
	0,                              /* tp_setattro */
	0,                              /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
	WithinImageBoundaryUP1D___doc__, /* tp_doc */
	0,                              /* tp_traverse */
	0,                              /* tp_clear */
	0,                              /* tp_richcompare */
	0,                              /* tp_weaklistoffset */
	0,                              /* tp_iter */
	0,                              /* tp_iternext */
	0,                              /* tp_methods */
	0,                              /* tp_members */
	0,                              /* tp_getset */
	&UnaryPredicate1D_Type,         /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	(initproc)WithinImageBoundaryUP1D___init__, /* tp_init */
	0,                              /* tp_alloc */
	0,                              /* tp_new */
};

///////////////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

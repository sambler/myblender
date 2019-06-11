/*
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
 */

/** \file
 * \ingroup freestyle
 */

#include "BPy_BlenderTextureShader.h"

#include "../../stroke/BasicStrokeShaders.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "../../../../python/generic/py_capi_utils.h"

///////////////////////////////////////////////////////////////////////////////////////////

//------------------------INSTANCE METHODS ----------------------------------

static char BlenderTextureShader___doc__[] =
    "Class hierarchy: :class:`freestyle.types.StrokeShader` > :class:`BlenderTextureShader`\n"
    "\n"
    "[Texture shader]\n"
    "\n"
    ".. method:: __init__(texture)\n"
    "\n"
    "   Builds a BlenderTextureShader object.\n"
    "\n"
    "   :arg texture: A line style texture slot or a shader node tree to define\n"
    "       a set of textures.\n"
    "   :type texture: :class:`bpy.types.LineStyleTextureSlot` or\n"
    "       :class:`bpy.types.ShaderNodeTree`\n"
    "\n"
    ".. method:: shade(stroke)\n"
    "\n"
    "   Assigns a blender texture slot to the stroke  shading in order to\n"
    "   simulate marks.\n"
    "\n"
    "   :arg stroke: A Stroke object.\n"
    "   :type stroke: :class:`freestyle.types.Stroke`\n";

static int BlenderTextureShader___init__(BPy_BlenderTextureShader *self,
                                         PyObject *args,
                                         PyObject *kwds)
{
  static const char *kwlist[] = {"texture", NULL};
  PyObject *obj;
  MTex *_mtex;
  bNodeTree *_nodetree;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", (char **)kwlist, &obj)) {
    return -1;
  }
  _mtex = (MTex *)PyC_RNA_AsPointer(obj, "LineStyleTextureSlot");
  if (_mtex) {
    self->py_ss.ss = new StrokeShaders::BlenderTextureShader(_mtex);
    return 0;
  }
  PyErr_Clear();
  _nodetree = (bNodeTree *)PyC_RNA_AsPointer(obj, "ShaderNodeTree");
  if (_nodetree) {
    self->py_ss.ss = new StrokeShaders::BlenderTextureShader(_nodetree);
    return 0;
  }
  PyErr_Format(PyExc_TypeError,
               "expected either 'LineStyleTextureSlot' or 'ShaderNodeTree', "
               "found '%.200s' instead",
               Py_TYPE(obj)->tp_name);
  return -1;
}

/*-----------------------BPy_BlenderTextureShader type definition ------------------------------*/

PyTypeObject BlenderTextureShader_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "BlenderTextureShader", /* tp_name */
    sizeof(BPy_BlenderTextureShader),                      /* tp_basicsize */
    0,                                                     /* tp_itemsize */
    0,                                                     /* tp_dealloc */
    0,                                                     /* tp_print */
    0,                                                     /* tp_getattr */
    0,                                                     /* tp_setattr */
    0,                                                     /* tp_reserved */
    0,                                                     /* tp_repr */
    0,                                                     /* tp_as_number */
    0,                                                     /* tp_as_sequence */
    0,                                                     /* tp_as_mapping */
    0,                                                     /* tp_hash  */
    0,                                                     /* tp_call */
    0,                                                     /* tp_str */
    0,                                                     /* tp_getattro */
    0,                                                     /* tp_setattro */
    0,                                                     /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,              /* tp_flags */
    BlenderTextureShader___doc__,                          /* tp_doc */
    0,                                                     /* tp_traverse */
    0,                                                     /* tp_clear */
    0,                                                     /* tp_richcompare */
    0,                                                     /* tp_weaklistoffset */
    0,                                                     /* tp_iter */
    0,                                                     /* tp_iternext */
    0,                                                     /* tp_methods */
    0,                                                     /* tp_members */
    0,                                                     /* tp_getset */
    &StrokeShader_Type,                                    /* tp_base */
    0,                                                     /* tp_dict */
    0,                                                     /* tp_descr_get */
    0,                                                     /* tp_descr_set */
    0,                                                     /* tp_dictoffset */
    (initproc)BlenderTextureShader___init__,               /* tp_init */
    0,                                                     /* tp_alloc */
    0,                                                     /* tp_new */
};

///////////////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

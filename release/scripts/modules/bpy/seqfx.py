# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8-80 compliant>

# for slightly faster access
from _bpy import seqfx as seqfx_module

# add = seqfx_module.add
#init = seqfx_module.dir
#poll = seqfx_module.poll
#call = seqfx_module.call
#as_string = seqfx_module.as_string
#get_rna = seqfx_module.get_rna
store_icu_range = seqfx_module.store_icu_range
get_default_fac = seqfx_module.get_default_fac

class BPySeqfx(object):
    '''
    Fake module like class.

     bpy.seqfx
    '''

    def __getattr__(self, module):
        '''
        gets a bpy.seqfx submodule
        '''
        if module.startswith('__'):
            raise AttributeError(module)
        return BPySeqfxSubMod(module)

    def __dir__(self):

        submodules = set()

        # add this classes functions
        for id_name in dir(self.__class__):
            if not id_name.startswith('__'):
                submodules.add(id_name)

        for id_name in dir():
            id_split = id_name.split('_OT_', 1)

            if len(id_split) == 2:
                submodules.add(id_split[0].lower())
            else:
                submodules.add(id_split[0])

        return list(submodules)

    def __repr__(self):
        return "<module like class 'bpy.seqfx'>"


class BPySeqfxSubMod(object):
    '''
    Utility class to fake submodules.

    eg. bpy.seqfx.object
    '''
    __keys__ = ('module',)

    def __init__(self, module):
        self.module = module

    def __getattr__(self, func):
        '''
        gets a bpy.seqfx.submodule function
        '''
        if func.startswith('__'):
            raise AttributeError(func)
        return BPySeqfxSubModOp(self.module, func)

    def __dir__(self):

        functions = set()

        module_upper = self.module.upper()

        for id_name in dir():
            id_split = id_name.split('_OT_', 1)
            if len(id_split) == 2 and module_upper == id_split[0]:
                functions.add(id_split[1])

        return list(functions)

    def __repr__(self):
        return "<module like class 'bpy.seqfx.%s'>" % self.module


class BPySeqfxSubModOp(object):
    '''
    Utility class to fake submodule operators.

    eg. bpy.seqfx.object.somefunc
    '''

    __keys__ = ('module', 'func')

    def _get_doc(self):
        return as_string(self.idname())

    @staticmethod
    def _parse_args(args):
        C_dict = None
        C_exec = 'EXEC_DEFAULT'

        if len(args) == 0:
            pass
        elif len(args) == 1:
            if type(args[0]) != str:
                C_dict = args[0]
            else:
                C_exec = args[0]
        elif len(args) == 2:
            C_exec, C_dict = args
        else:
            raise ValueError("1 or 2 args execution context is supported")

        return C_dict, C_exec

    @staticmethod
    def _scene_update(context):
        scene = context.scene
        if scene:  # None in backgroud mode
            scene.update()
        else:
            import bpy
            for scene in bpy.data.scenes:
                scene.update()

    __doc__ = property(_get_doc)

    def __init__(self, module, func):
        self.module = module
        self.func = func

    def poll(self, *args):
        C_dict, C_exec = BPySeqfxSubModOp._parse_args(args)
        return poll(self.idname_py(), C_dict, C_exec)

    def idname(self):
        # submod.foo -> SUBMOD_OT_foo
        return self.module.upper() + "_OT_" + self.func

    def idname_py(self):
        # submod.foo -> SUBMOD_OT_foo
        return self.module + "." + self.func

    def __call__(self, *args, **kw):
        import bpy
        context = bpy.context

        # Get the operator from blender
        wm = context.window_manager

        # run to account for any rna values the user changes.
        BPySeqfxSubModOp._scene_update(context)

        if args:
            C_dict, C_exec = BPySeqfxSubModOp._parse_args(args)
            ret = call(self.idname_py(), C_dict, kw, C_exec)
        else:
            ret = call(self.idname_py(), None, kw)

        if 'FINISHED' in ret and context.window_manager == wm:
            BPySeqfxSubModOp._scene_update(context)

        return ret

    def get_rna(self):
        '''
        currently only used for 'bl_rna'
        '''
        return get_rna(self.idname())

    def __repr__(self):  # useful display, repr(op)
        import bpy
        idname = self.idname()
        as_string = as_string(idname)
        seqfx_class = getattr(bpy.types, idname)
        descr = seqfx_class.bl_rna.description
        # XXX, workaround for not registering
        # every __doc__ to save time on load.
        if not descr:
            descr = seqfx_class.__doc__
            if not descr:
                descr = ""

        return "# %s\n%s" % (descr, as_string)

    def __str__(self):  # used for print(...)
        return "<function bpy.seqfx.%s.%s at 0x%x'>" % \
                (self.module, self.func, id(self))


seqfx_fake_module = BPySeqfx()

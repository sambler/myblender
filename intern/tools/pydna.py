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

# <pep8 compliant>

"""
Experemental module (for developers only), exposes DNA data via python
uses no blender/python modules, pure python + autogenerated ctypes api.

Exposes:
 * pydna.main: main database access.
 * pydna.types: a module of all DNA ctypes structures.

 * Utility method: CAST(dna_type)
 * Utility method for ListBase: ITER(dna_type)

Example:
 import sys
 sys.path.append("/b/intern/tools")

 import pydna

 for obj in pydna.main.object.ITER("Object"):
     print("Object:", obj.id.name)
     if obj.data:
         data = pydna.types.ID.from_address(obj.data)
         print("  ObData:", data.name)
"""

import ctypes
import struct


def _api():

    def is_ctypes_subclass(other, main_class):
        while other:
            if other is main_class:
                return True
            other = getattr(other, "_type_", None)
        return False

    def is_ctypes_initialized(other):
        _other = other
        while other:
            if hasattr(other, "_fields_"):
                return True
            other = getattr(other, "_type_", None)
        print(_other, "NOT INIT")
        return False

    def is_ctypes_base(other):
        while type(getattr(other, "_type_", "")) != str:
            other = other._type_
        return other
    
    class MixIn:
        pass
    
    blend_cdll = ctypes.CDLL("")
    blend_lib = ctypes.LibraryLoader("")
    
    def blend_parse_dna():
        # from dna.c
        sdna_str_pt = blend_cdll.DNAstr
        sdna_len_pt = blend_cdll.DNAlen

        # cast
        sdna_len_pt = ctypes.c_void_p.from_address(ctypes.addressof(sdna_len_pt))
        sdna_len = ctypes.c_int.from_address(sdna_len_pt.value)

        blend_sdna = ctypes.string_at(sdna_str_pt, sdna_len)
        
        ofs = 0
        assert(blend_sdna[ofs:ofs + 8] == b'SDNANAME')
        ofs += 8

        sdna_names_len = struct.unpack("i", blend_sdna[ofs:ofs + 4])[0]
        ofs += 4
        
        blend_sdna_names = blend_sdna[ofs:].split(b'\0', sdna_names_len)
        blend_sdna_remainder = blend_sdna_names.pop(-1)  # last item is not a name.
        ofs = len(blend_sdna) - len(blend_sdna_remainder)
        ofs = (ofs + 3) & ~3

        assert(blend_sdna[ofs:ofs + 4] == b'TYPE')
        ofs += 4

        sdna_types_len = struct.unpack("i", blend_sdna[ofs:ofs + 4])[0]
        ofs += 4
        
        blend_sdna_types = blend_sdna[ofs:].split(b'\0', sdna_types_len)
        blend_sdna_remainder = blend_sdna_types.pop(-1)
        ofs = len(blend_sdna) - len(blend_sdna_remainder)
        ofs = (ofs + 3) & ~3

        assert(blend_sdna[ofs:ofs + 4] == b'TLEN')
        ofs += 4

        blend_sdna_typelens = struct.unpack("%dh" % sdna_types_len, blend_sdna[ofs:ofs + (sdna_types_len * 2)])
        ofs += sdna_types_len * 2
        ofs = (ofs + 3) & ~3

        # array of pointers to short arrays
        assert(blend_sdna[ofs:ofs + 4] == b'STRC')
        ofs += 4
        
        sdna_structs_len = struct.unpack("i", blend_sdna[ofs:ofs + 4])[0]
        ofs += 4

        blend_sdna_structs = []

        for i in range(sdna_structs_len):
            struct_type, struct_tot = struct.unpack("hh", blend_sdna[ofs:ofs + 4])
            ofs += 4
            struct_type_name_pairs = struct.unpack("%dh" % struct_tot * 2, blend_sdna[ofs:ofs + (struct_tot * 4)])
        
            # convert into pairs, easier to understand (type, name)
            struct_type_name_pairs = [(struct_type_name_pairs[j], struct_type_name_pairs[j + 1]) for j in range(0, struct_tot * 2, 2)]

            blend_sdna_structs.append((struct_type, struct_type_name_pairs))
            ofs += struct_tot * 4
            # ofs = (ofs + 1) & ~1
        
        return blend_sdna_names, blend_sdna_types, blend_sdna_typelens, blend_sdna_structs

    def create_dna_structs(blend_sdna_names, blend_sdna_types, blend_sdna_typelens, blend_sdna_structs):
        
        # create all subclasses of ctypes.Structure
        ctypes_structs = {name: type(name.decode(), (ctypes.Structure, MixIn), {}) for name in blend_sdna_types}
        ctypes_basic = {b"float": ctypes.c_float, b"double": ctypes.c_double, b"int": ctypes.c_int, b"short": ctypes.c_short, b"char": ctypes.c_char, b"void": ctypes.c_void_p}
        ctypes_fields = {}
        
        # collect fields
        for struct_id, struct_type_name_pairs in blend_sdna_structs:
            struct_name = blend_sdna_types[struct_id]
            ctype_struct = ctypes_structs[struct_name]
            fields = []
            
            for stype, sname in struct_type_name_pairs:
                name_string = blend_sdna_names[sname]
                type_string = blend_sdna_types[stype]
                type_py = ctypes_basic.get(type_string)
                if type_py is None:
                    type_py = ctypes_structs.get(type_string)

                # todo, these might need to be changed
                name_string = name_string.replace(b"(", b"")
                name_string = name_string.replace(b")", b"")

                # * First parse the pointer *
                pointer_count = 0
                while name_string[0] == 42:  # '*'
                    pointer_count += 1
                    name_string = name_string[1:]
                
                # alredy a pointer
                if type_py is ctypes.c_void_p:
                    pointer_count -= 1
                elif type_py is ctypes.c_char and pointer_count == 1:
                    type_py = ctypes.c_char_p
                    pointer_count = 0
                
                if pointer_count < 0:
                    Exception("error parsing pointer")

                for i in range(pointer_count):
                    type_py = ctypes.POINTER(type_py)
                
                # * Now parse the array [] *
                if b'[' in name_string:
                    name_string = name_string.replace(b'[', b' ')
                    name_string = name_string.replace(b']', b' ')
                    name_split = name_string.split()
                    name_string = name_split[0]
                    for array_dim in reversed(name_split[1:]):
                        type_py = type_py * int(array_dim)

                fields.append((name_string.decode(), type_py))

            ctypes_fields[struct_name] = fields
        
        # apply fields all in one go!
        for struct_id, struct_type_name_pairs in blend_sdna_structs:
            struct_name = blend_sdna_types[struct_id]
            ctype_struct = ctypes_structs[struct_name]
            try:
                ctype_struct._fields_ = ctypes_fields[struct_name]
            except:
                print("Error:", struct_name)
                import traceback
                traceback.print_exc()
                # print(fields)

        # test fields
        for struct_id, struct_type_name_pairs in blend_sdna_structs:
            ctype_struct = ctypes_structs[blend_sdna_types[struct_id]]
            if blend_sdna_typelens[struct_id] != ctypes.sizeof(ctype_struct):
                print("Size Mismatch for %r blender:%d vs python:%d" % (blend_sdna_types[struct_id], blend_sdna_typelens[struct_id], ctypes.sizeof(ctype_struct)))

        return ctypes_structs

    def decorate_api(struct_dict):
        
        # * Decotate the api *
        
        # listbase iter
        
        type_cast_lb = struct_dict[b'ListBase']
        type_cast_link = struct_dict[b'Link']

        def list_base_iter(self, type_name):
            type_cast = struct_dict[type_name.encode('ASCII')]
            try:
                ret = type_cast_link.from_address(ctypes.addressof(self.first))
            except:
                ret = type_cast_link.from_address(self.first)

            while ret is not None:
                return_value = type_cast.from_address(ctypes.addressof(ret))
                try:
                    next_pointer = getattr(ret.next, "contents")
                except:
                    next_pointer = None

                if next_pointer:
                    ret = type_cast_link.from_address(ctypes.addressof(next_pointer))
                else:
                    ret = None

                yield return_value

        struct_dict[b'ListBase'].ITER = list_base_iter
        
        def CAST(self, to):
            type_cast = struct_dict[to.encode('ASCII')]
            return type_cast.from_address(ctypes.addressof(self))
        
        MixIn.CAST = CAST

    blend_sdna_names, blend_sdna_types, blend_sdna_typelens, blend_sdna_structs = blend_parse_dna()
    struct_dict = create_dna_structs(blend_sdna_names, blend_sdna_types, blend_sdna_typelens, blend_sdna_structs)

    # print out all structs
    '''
    for struct_id, struct_type_name_pairs in blend_sdna_structs:
        print("")
        sruct_name = blend_sdna_types[struct_id].decode()
        print("typedef struct %s {" % sruct_name)
        for stype, sname in struct_type_name_pairs:
            print("    %s %s;" % (blend_sdna_types[stype].decode(), blend_sdna_names[sname].decode()))
        print("} %s;" % sruct_name)
    '''
    
    decorate_api(struct_dict)  # not essential but useful
    
    # manually wrap Main
    Main = type("Main", (ctypes.Structure, ), {})
    _lb = struct_dict[b"ListBase"]
    Main._fields_ = [("next", ctypes.POINTER(Main)),
                     ("prev", ctypes.POINTER(Main)),
                     ("name", ctypes.c_char * 240),
                     ("versionfile", ctypes.c_short),
                     ("subversionfile", ctypes.c_short),
                     ("minversionfile", ctypes.c_short),
                     ("minsubversionfile", ctypes.c_short),
                     ("curlib", ctypes.POINTER(struct_dict[b"Library"])),
                     ("scene", _lb),
                     ("library", _lb),
                     ("object", _lb),
                     ("mesh", _lb),
                     ("curve", _lb),
                     ("mball", _lb),
                     ("mat", _lb),
                     ("tex", _lb),
                     ("image", _lb),
                     ("latt", _lb),
                     ("lamp", _lb),
                     ("camera", _lb),
                     ("ipo", _lb),
                     ("key", _lb),
                     ("world", _lb),
                     ("screen", _lb),
                     ("script", _lb),
                     ("vfont", _lb),
                     ("text", _lb),
                     ("sound", _lb),
                     ("group", _lb),
                     ("armature", _lb),
                     ("action", _lb),
                     ("nodetree", _lb),
                     ("brush", _lb),
                     ("particle", _lb),
                     ("wm", _lb),
                     ("gpencil", _lb),
                     ]
    del _lb

    # import bpy
    # main = Main.from_address(bpy.data.as_pointer())
    # main is the first pointer in Global.
    main_address = ctypes.POINTER(ctypes.c_void_p).from_address(ctypes.addressof(blend_cdll.G)).contents.value
    main = Main.from_address(main_address)
    
    return main, struct_dict

main, _struct_dict = _api()

# types dict
types = type(ctypes)("pydna.types")
types.__dict__.update({s.__name__: s for s in _struct_dict.values()})

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
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 */

/** \file
 * \ingroup DNA
 * \brief ID and Library types, which are fundamental for sdna.
 */

#ifndef __DNA_ID_H__
#define __DNA_ID_H__

#include "DNA_listBase.h"

#ifdef __cplusplus
extern "C" {
#endif

struct FileData;
struct GPUTexture;
struct ID;
struct Library;
struct PackedFile;

/* Runtime display data */
struct DrawData;
typedef void (*DrawDataInitCb)(struct DrawData *engine_data);
typedef void (*DrawDataFreeCb)(struct DrawData *engine_data);

#
#
typedef struct DrawData {
  struct DrawData *next, *prev;
  struct DrawEngineType *engine_type;
  /* Only nested data, NOT the engine data itself. */
  DrawDataFreeCb free;
  /* Accumulated recalc flags, which corresponds to ID->recalc flags. */
  int recalc;
} DrawData;

typedef struct DrawDataList {
  struct DrawData *first, *last;
} DrawDataList;

typedef struct IDPropertyData {
  void *pointer;
  ListBase group;
  /** Note, we actually fit a double into these two ints. */
  int val, val2;
} IDPropertyData;

typedef struct IDProperty {
  struct IDProperty *next, *prev;
  char type, subtype;
  short flag;
  /** MAX_IDPROP_NAME. */
  char name[64];

  /* saved is used to indicate if this struct has been saved yet.
   * seemed like a good idea as a '_pad' var was needed anyway :) */
  int saved;
  /** Note, alignment for 64 bits. */
  IDPropertyData data;

  /* array length, also (this is important!) string length + 1.
   * the idea is to be able to reuse array realloc functions on strings.*/
  int len;

  /* Strings and arrays are both buffered, though the buffer isn't saved. */
  /* totallen is total length of allocated array/string, including a buffer.
   * Note that the buffering is mild; the code comes from python's list implementation. */
  int totallen;
} IDProperty;

#define MAX_IDPROP_NAME 64
#define DEFAULT_ALLOC_FOR_NULL_STRINGS 64

/*->type*/
enum {
  IDP_STRING = 0,
  IDP_INT = 1,
  IDP_FLOAT = 2,
  IDP_ARRAY = 5,
  IDP_GROUP = 6,
  IDP_ID = 7,
  IDP_DOUBLE = 8,
  IDP_IDPARRAY = 9,
  IDP_NUMTYPES = 10,
};

/*->subtype */

/* IDP_STRING */
enum {
  IDP_STRING_SUB_UTF8 = 0, /* default */
  IDP_STRING_SUB_BYTE = 1, /* arbitrary byte array, _not_ null terminated */
};

/* IDP_GROUP */
enum {
  /** Default. */
  IDP_GROUP_SUB_NONE = 0,
  /** Object mode settings. */
  IDP_GROUP_SUB_MODE_OBJECT = 1,
  /** Mesh edit mode settings. */
  IDP_GROUP_SUB_MODE_EDIT = 2,
  /** Render engine settings. */
  IDP_GROUP_SUB_ENGINE_RENDER = 3,
  /** Data override. */
  IDP_GROUP_SUB_OVERRIDE = 4,
  /** Weight paint mode settings. */
  IDP_GROUP_SUB_MODE_PAINT_WEIGHT = 5,
  /** Vertex paint mode settings. */
  IDP_GROUP_SUB_MODE_PAINT_VERTEX = 6,
};

/*->flag*/
enum {
  /** This IDProp may be statically overridden.
   * Should only be used/be relevant for custom properties. */
  IDP_FLAG_OVERRIDABLE_STATIC = 1 << 0,

  /** This means the property is set but RNA will return false when checking
   * 'RNA_property_is_set', currently this is a runtime flag */
  IDP_FLAG_GHOST = 1 << 7,
};

/* add any future new id property types here.*/

/* Static ID override structs. */

typedef struct IDOverrideStaticPropertyOperation {
  struct IDOverrideStaticPropertyOperation *next, *prev;

  /* Type of override. */
  short operation;
  short flag;
  char _pad0[4];

  /* Sub-item references, if needed (for arrays or collections only).
   * We need both reference and local values to allow e.g. insertion into collections
   * (constraints, modifiers...).
   * In collection case, if names are defined, they are used in priority.
   * Names are pointers (instead of char[64]) to save some space, NULL when unset.
   * Indices are -1 when unset. */
  char *subitem_reference_name;
  char *subitem_local_name;
  int subitem_reference_index;
  int subitem_local_index;
} IDOverrideStaticPropertyOperation;

/* IDOverridePropertyOperation->operation. */
enum {
  /* Basic operations. */
  IDOVERRIDESTATIC_OP_NOOP = 0, /* Special value, forbids any overriding. */

  IDOVERRIDESTATIC_OP_REPLACE = 1, /* Fully replace local value by reference one. */

  /* Numeric-only operations. */
  IDOVERRIDESTATIC_OP_ADD = 101, /* Add local value to reference one. */
  /* Subtract local value from reference one (needed due to unsigned values etc.). */
  IDOVERRIDESTATIC_OP_SUBTRACT = 102,
  /* Multiply reference value by local one (more useful than diff for scales and the like). */
  IDOVERRIDESTATIC_OP_MULTIPLY = 103,

  /* Collection-only operations. */
  IDOVERRIDESTATIC_OP_INSERT_AFTER = 201,  /* Insert after given reference's subitem. */
  IDOVERRIDESTATIC_OP_INSERT_BEFORE = 202, /* Insert before given reference's subitem. */
  /* We can add more if needed (move, delete, ...). */
};

/* IDOverridePropertyOperation->flag. */
enum {
  /** User cannot remove that override operation. */
  IDOVERRIDESTATIC_FLAG_MANDATORY = 1 << 0,
  /** User cannot change that override operation. */
  IDOVERRIDESTATIC_FLAG_LOCKED = 1 << 1,
};

/** A single overridden property, contain all operations on this one. */
typedef struct IDOverrideStaticProperty {
  struct IDOverrideStaticProperty *next, *prev;

  /**
   * Path from ID to overridden property.
   * *Does not* include indices/names for final arrays/collections items.
   */
  char *rna_path;

  /** List of overriding operations (IDOverridePropertyOperation) applied to this property. */
  ListBase operations;
} IDOverrideStaticProperty;

/* Main container for all overriding data info of a data-block. */
typedef struct IDOverrideStatic {
  /** Reference linked ID which this one overrides. */
  struct ID *reference;
  /** List of IDOverrideProperty structs. */
  ListBase properties;

  short flag;
  char _pad[6];

  /* Read/write data. */
  /* Temp ID storing extra override data (used for differential operations only currently).
   * Always NULL outside of read/write context. */
  struct ID *storage;
} IDOverrideStatic;

enum eStaticOverride_Flag {
  STATICOVERRIDE_AUTO = 1 << 0, /* Allow automatic generation of overriding rules. */
};

/* watch it: Sequence has identical beginning. */
/**
 * ID is the first thing included in all serializable types. It
 * provides a common handle to place all data in double-linked lists.
 * */

/* 2 characters for ID code and 64 for actual name */
#define MAX_ID_NAME 66

/* There's a nasty circular dependency here.... 'void *' to the rescue! I
 * really wonder why this is needed. */
typedef struct ID {
  void *next, *prev;
  struct ID *newid;
  struct Library *lib;
  /** MAX_ID_NAME. */
  char name[66];
  /**
   * LIB_... flags report on status of the datablock this ID belongs to
   * (persistent, saved to and read from .blend).
   */
  short flag;
  /**
   * LIB_TAG_... tags (runtime only, cleared at read time).
   */
  int tag;
  int us;
  int icon_id;
  int recalc;
  char _pad[4];
  IDProperty *properties;

  /** Reference linked ID which this one overrides. */
  IDOverrideStatic *override_static;

  /**
   * Only set for datablocks which are coming from copy-on-write, points to
   * the original version of it.
   */
  struct ID *orig_id;

  void *py_instance;
} ID;

/**
 * For each library file used, a Library struct is added to Main
 * WARNING: readfile.c, expand_doit() reads this struct without DNA check!
 */
typedef struct Library {
  ID id;
  struct FileData *filedata;
  /** Path name used for reading, can be relative and edited in the outliner. */
  char name[1024];

  /**
   * Absolute filepath, this is only for convenience,
   * 'name' is the real path used on file read but in
   * some cases its useful to access the absolute one.
   * This is set on file read.
   * Use BKE_library_filepath_set() rather than setting 'name'
   * directly and it will be kept in sync - campbell */
  char filepath[1024];

  /** Set for indirectly linked libs, used in the outliner and while reading. */
  struct Library *parent;

  struct PackedFile *packedfile;

  /* Temp data needed by read/write code. */
  int temp_index;
  /** See BLENDER_VERSION, BLENDER_SUBVERSION, needed for do_versions. */
  short versionfile, subversionfile;
} Library;

enum eIconSizes {
  ICON_SIZE_ICON = 0,
  ICON_SIZE_PREVIEW = 1,

  NUM_ICON_SIZES,
};

/* for PreviewImage->flag */
enum ePreviewImage_Flag {
  PRV_CHANGED = (1 << 0),
  PRV_USER_EDITED = (1 << 1), /* if user-edited, do not auto-update this anymore! */
};

/* for PreviewImage->tag */
enum {
  PRV_TAG_DEFFERED = (1 << 0),           /* Actual loading of preview is deferred. */
  PRV_TAG_DEFFERED_RENDERING = (1 << 1), /* Deferred preview is being loaded. */
  PRV_TAG_DEFFERED_DELETE = (1 << 2),    /* Deferred preview should be deleted asap. */
};

typedef struct PreviewImage {
  /* All values of 2 are really NUM_ICON_SIZES */
  unsigned int w[2];
  unsigned int h[2];
  short flag[2];
  short changed_timestamp[2];
  unsigned int *rect[2];

  /* Runtime-only data. */
  struct GPUTexture *gputexture[2];
  /** Used by previews outside of ID context. */
  int icon_id;

  /** Runtime data. */
  short tag;
  char _pad[2];
} PreviewImage;

#define PRV_DEFERRED_DATA(prv) \
  (CHECK_TYPE_INLINE(prv, PreviewImage *), \
   BLI_assert((prv)->tag & PRV_TAG_DEFFERED), \
   (void *)((prv) + 1))

/**
 * Defines for working with IDs.
 *
 * The tags represent types! This is a dirty way of enabling RTTI. The
 * sig_byte end endian defines aren't really used much.
 */

#ifdef __BIG_ENDIAN__
/* big endian */
#  define MAKE_ID2(c, d) ((c) << 8 | (d))
#else
/* little endian  */
#  define MAKE_ID2(c, d) ((d) << 8 | (c))
#endif

/**
 * ID from database.
 *
 * Written to #BHead.code (for file IO)
 * and the first 2 bytes of #ID.name (for runtime checks, see #GS macro).
 */
typedef enum ID_Type {
  ID_SCE = MAKE_ID2('S', 'C'), /* Scene */
  ID_LI = MAKE_ID2('L', 'I'),  /* Library */
  ID_OB = MAKE_ID2('O', 'B'),  /* Object */
  ID_ME = MAKE_ID2('M', 'E'),  /* Mesh */
  ID_CU = MAKE_ID2('C', 'U'),  /* Curve */
  ID_MB = MAKE_ID2('M', 'B'),  /* MetaBall */
  ID_MA = MAKE_ID2('M', 'A'),  /* Material */
  ID_TE = MAKE_ID2('T', 'E'),  /* Tex (Texture) */
  ID_IM = MAKE_ID2('I', 'M'),  /* Image */
  ID_LT = MAKE_ID2('L', 'T'),  /* Lattice */
  ID_LA = MAKE_ID2('L', 'A'),  /* Light */
  ID_CA = MAKE_ID2('C', 'A'),  /* Camera */
  ID_IP = MAKE_ID2('I', 'P'),  /* Ipo (depreciated, replaced by FCurves) */
  ID_KE = MAKE_ID2('K', 'E'),  /* Key (shape key) */
  ID_WO = MAKE_ID2('W', 'O'),  /* World */
  ID_SCR = MAKE_ID2('S', 'R'), /* Screen */
  ID_VF = MAKE_ID2('V', 'F'),  /* VFont (Vector Font) */
  ID_TXT = MAKE_ID2('T', 'X'), /* Text */
  ID_SPK = MAKE_ID2('S', 'K'), /* Speaker */
  ID_SO = MAKE_ID2('S', 'O'),  /* Sound */
  ID_GR = MAKE_ID2('G', 'R'),  /* Group */
  ID_AR = MAKE_ID2('A', 'R'),  /* bArmature */
  ID_AC = MAKE_ID2('A', 'C'),  /* bAction */
  ID_NT = MAKE_ID2('N', 'T'),  /* bNodeTree */
  ID_BR = MAKE_ID2('B', 'R'),  /* Brush */
  ID_PA = MAKE_ID2('P', 'A'),  /* ParticleSettings */
  ID_GD = MAKE_ID2('G', 'D'),  /* bGPdata, (Grease Pencil) */
  ID_WM = MAKE_ID2('W', 'M'),  /* WindowManager */
  ID_MC = MAKE_ID2('M', 'C'),  /* MovieClip */
  ID_MSK = MAKE_ID2('M', 'S'), /* Mask */
  ID_LS = MAKE_ID2('L', 'S'),  /* FreestyleLineStyle */
  ID_PAL = MAKE_ID2('P', 'L'), /* Palette */
  ID_PC = MAKE_ID2('P', 'C'),  /* PaintCurve  */
  ID_CF = MAKE_ID2('C', 'F'),  /* CacheFile */
  ID_WS = MAKE_ID2('W', 'S'),  /* WorkSpace */
  ID_LP = MAKE_ID2('L', 'P'),  /* LightProbe */
} ID_Type;

/* Only used as 'placeholder' in .blend files for directly linked datablocks. */
#define ID_LINK_PLACEHOLDER MAKE_ID2('I', 'D') /* (internal use only) */

/* Deprecated. */
#define ID_SCRN MAKE_ID2('S', 'N')

/* NOTE! Fake IDs, needed for g.sipo->blocktype or outliner */
#define ID_SEQ MAKE_ID2('S', 'Q')
/* constraint */
#define ID_CO MAKE_ID2('C', 'O')
/* pose (action channel, used to be ID_AC in code, so we keep code for backwards compat) */
#define ID_PO MAKE_ID2('A', 'C')
/* used in outliner... */
#define ID_NLA MAKE_ID2('N', 'L')
/* fluidsim Ipo */
#define ID_FLUIDSIM MAKE_ID2('F', 'S')

#define ID_FAKE_USERS(id) ((((ID *)id)->flag & LIB_FAKEUSER) ? 1 : 0)
#define ID_REAL_USERS(id) (((ID *)id)->us - ID_FAKE_USERS(id))
#define ID_EXTRA_USERS(id) (((ID *)id)->tag & LIB_TAG_EXTRAUSER ? 1 : 0)

#define ID_CHECK_UNDO(id) \
  ((GS((id)->name) != ID_SCR) && (GS((id)->name) != ID_WM) && (GS((id)->name) != ID_WS))

#define ID_BLEND_PATH(_bmain, _id) \
  ((_id)->lib ? (_id)->lib->filepath : BKE_main_blendfile_path((_bmain)))
#define ID_BLEND_PATH_FROM_GLOBAL(_id) \
  ((_id)->lib ? (_id)->lib->filepath : BKE_main_blendfile_path_from_global())

#define ID_MISSING(_id) (((_id)->tag & LIB_TAG_MISSING) != 0)

#define ID_IS_LINKED(_id) (((ID *)(_id))->lib != NULL)

#define ID_IS_STATIC_OVERRIDE(_id) \
  (((ID *)(_id))->override_static != NULL && ((ID *)(_id))->override_static->reference != NULL)

#define ID_IS_STATIC_OVERRIDE_TEMPLATE(_id) \
  (((ID *)(_id))->override_static != NULL && ((ID *)(_id))->override_static->reference == NULL)

#define ID_IS_STATIC_OVERRIDE_AUTO(_id) \
  (!ID_IS_LINKED((_id)) && ID_IS_STATIC_OVERRIDE((_id)) && \
   (((ID *)(_id))->override_static->flag & STATICOVERRIDE_AUTO))

/* No copy-on-write for these types.
 * Keep in sync with check_datablocks_copy_on_writable and deg_copy_on_write_is_needed */
#define ID_TYPE_IS_COW(_id_type) (!ELEM(_id_type, ID_BR, ID_LS, ID_PAL, ID_IM))

#ifdef GS
#  undef GS
#endif
#define GS(a) \
  (CHECK_TYPE_ANY(a, char *, const char *, char[66], const char[66]), \
   (ID_Type)(*((const short *)(a))))

#define ID_NEW_SET(_id, _idn) \
  (((ID *)(_id))->newid = (ID *)(_idn), \
   ((ID *)(_id))->newid->tag |= LIB_TAG_NEW, \
   (void *)((ID *)(_id))->newid)
#define ID_NEW_REMAP(a) \
  if ((a) && (a)->id.newid) \
  (a) = (void *)(a)->id.newid

/* id->flag (persitent). */
enum {
  LIB_FAKEUSER = 1 << 9,
};

/**
 * id->tag (runtime-only).
 *
 * Those flags belong to three different categories,
 * which have different expected handling in code:
 *
 * - RESET_BEFORE_USE: piece of code that wants to use such flag
 *   has to ensure they are properly 'reset' first.
 * - RESET_AFTER_USE: piece of code that wants to use such flag has to ensure they are properly
 *   'reset' after usage
 *   (though 'lifetime' of those flags is a bit fuzzy, e.g. _RECALC ones are reset on depsgraph
 *   evaluation...).
 * - RESET_NEVER: those flags are 'status' one, and never actually need any reset
 *   (except on initialization during .blend file reading).
 */
enum {
  /* RESET_NEVER Datablock is from current .blend file. */
  LIB_TAG_LOCAL = 0,
  /* RESET_NEVER Datablock is from a library,
   * but is used (linked) directly by current .blend file. */
  LIB_TAG_EXTERN = 1 << 0,
  /* RESET_NEVER Datablock is from a library,
   * and is only used (linked) inderectly through other libraries. */
  LIB_TAG_INDIRECT = 1 << 1,

  /* RESET_AFTER_USE Flag used internally in readfile.c,
   * to mark IDs needing to be expanded (only done once). */
  LIB_TAG_NEED_EXPAND = 1 << 3,
  /* RESET_AFTER_USE Flag used internally in readfile.c to mark ID
   * placeholders for linked datablocks needing to be read. */
  LIB_TAG_ID_LINK_PLACEHOLDER = 1 << 4,
  /* RESET_AFTER_USE */
  LIB_TAG_NEED_LINK = 1 << 5,

  /* RESET_NEVER tag datablock as a place-holder
   * (because the real one could not be linked from its library e.g.). */
  LIB_TAG_MISSING = 1 << 6,

  /* RESET_NEVER tag datablock as being up-to-date regarding its reference. */
  LIB_TAG_OVERRIDESTATIC_REFOK = 1 << 9,
  /* RESET_NEVER tag datablock as needing an auto-override execution, if enabled. */
  LIB_TAG_OVERRIDESTATIC_AUTOREFRESH = 1 << 17,

  /* tag datablock has having an extra user. */
  LIB_TAG_EXTRAUSER = 1 << 2,
  /* tag datablock has having actually increased usercount for the extra virtual user. */
  LIB_TAG_EXTRAUSER_SET = 1 << 7,

  /* RESET_AFTER_USE tag newly duplicated/copied IDs.
   * Also used internally in readfile.c to mark datablocks needing do_versions. */
  LIB_TAG_NEW = 1 << 8,
  /* RESET_BEFORE_USE free test flag.
   * TODO make it a RESET_AFTER_USE too. */
  LIB_TAG_DOIT = 1 << 10,
  /* RESET_AFTER_USE tag existing data before linking so we know what is new. */
  LIB_TAG_PRE_EXISTING = 1 << 11,

  /* The datablock is a copy-on-write/localized version. */
  LIB_TAG_COPIED_ON_WRITE = 1 << 12,
  LIB_TAG_COPIED_ON_WRITE_EVAL_RESULT = 1 << 13,
  LIB_TAG_LOCALIZED = 1 << 14,

  /* RESET_NEVER tag datablock for freeing etc. behavior
   * (usually set when copying real one into temp/runtime one). */
  LIB_TAG_NO_MAIN = 1 << 15,          /* Datablock is not listed in Main database. */
  LIB_TAG_NO_USER_REFCOUNT = 1 << 16, /* Datablock does not refcount usages of other IDs. */
  /* Datablock was not allocated by standard system (BKE_libblock_alloc), do not free its memory
   * (usual type-specific freeing is called though). */
  LIB_TAG_NOT_ALLOCATED = 1 << 18,
};

/* Tag given ID for an update in all the dependency graphs. */
typedef enum IDRecalcFlag {
  /***************************************************************************
   * Individual update tags, this is what ID gets tagged for update with. */

  /* ** Object transformation changed. ** */
  ID_RECALC_TRANSFORM = (1 << 0),

  /* ** Object geometry changed. **
   *
   * When object of armature type gets tagged with this flag, it's pose is
   * re-evaluated.
   * When object of other type is tagged with this flag it makes the modifier
   * stack to be re-evaluated.
   * When object data type (mesh, curve, ...) gets tagged with this flag it
   * makes all objects which shares this datablock to be updated. */
  ID_RECALC_GEOMETRY = (1 << 1),

  /* ** Animation or time changed and animation is to be re-evaluated. ** */
  ID_RECALC_ANIMATION = (1 << 2),

  /* ** Particle system changed. ** */
  /* Only do pathcache etc. */
  ID_RECALC_PSYS_REDO = (1 << 3),
  /* Reset everything including pointcache. */
  ID_RECALC_PSYS_RESET = (1 << 4),
  /* Only child settings changed. */
  ID_RECALC_PSYS_CHILD = (1 << 5),
  /* Physics type changed. */
  ID_RECALC_PSYS_PHYS = (1 << 6),

  /* ** Material and shading ** */

  /* For materials and node trees this means that topology of the shader tree
   * changed, and the shader is to be recompiled.
   * For objects it means that the draw batch cache is to be redone. */
  ID_RECALC_SHADING = (1 << 7),
  /* TODO(sergey): Consider adding an explicit ID_RECALC_SHADING_PARAMATERS
   * which can be used for cases when only socket value changed, to speed up
   * redraw update in that case. */

  /* Selection of the ID itself or its components (for example, vertices) did
   * change, and all the drawing data is to eb updated. */
  ID_RECALC_SELECT = (1 << 9),
  /* Flags on the base did change, and is to be compied onto all the copies of
   * corresponding objects. */
  ID_RECALC_BASE_FLAGS = (1 << 10),
  ID_RECALC_POINT_CACHE = (1 << 11),
  /* Only inform editors about the change. Is used to force update of editors
   * when datablock which is not a part of dependency graph did change.
   *
   * For example, brush texture did change and the preview is to be
   * re-rendered. */
  ID_RECALC_EDITORS = (1 << 12),

  /* ** Update copy on write component. **
   * This is most generic tag which should only be used when nothing else
   * matches.
   */
  ID_RECALC_COPY_ON_WRITE = (1 << 13),

  /* Sequences in the sequencer did change.
   * Use this tag with a scene ID which owns the sequences. */
  ID_RECALC_SEQUENCER_STRIPS = (1 << 14),

  ID_RECALC_AUDIO_SEEK = (1 << 15),
  ID_RECALC_AUDIO_FPS = (1 << 16),
  ID_RECALC_AUDIO_VOLUME = (1 << 17),
  ID_RECALC_AUDIO_MUTE = (1 << 18),
  ID_RECALC_AUDIO_LISTENER = (1 << 19),

  /***************************************************************************
   * Pseudonyms, to have more semantic meaning in the actual code without
   * using too much low-level and implementation specific tags. */

  /* Update animation datablock itself, without doing full re-evaluation of
   * all dependent objects. */
  ID_RECALC_ANIMATION_NO_FLUSH = ID_RECALC_COPY_ON_WRITE,

  /***************************************************************************
   * Aggregate flags, use only for checks on runtime.
   * Do NOT use those for tagging. */

  /* Identifies that SOMETHING has been changed in this ID. */
  ID_RECALC_ALL = ~(0),
  /* Identifies that something in particle system did change. */
  ID_RECALC_PSYS_ALL = (ID_RECALC_PSYS_REDO | ID_RECALC_PSYS_RESET | ID_RECALC_PSYS_CHILD |
                        ID_RECALC_PSYS_PHYS),

} IDRecalcFlag;

/* To filter ID types (filter_id) */
/* XXX We cannot put all needed IDs inside an enum...
 *     We'll have to see whether we can fit all needed ones inside 32 values,
 *     or if we need to fallback to longlong defines :/
 */
enum {
  FILTER_ID_AC = (1 << 0),
  FILTER_ID_AR = (1 << 1),
  FILTER_ID_BR = (1 << 2),
  FILTER_ID_CA = (1 << 3),
  FILTER_ID_CU = (1 << 4),
  FILTER_ID_GD = (1 << 5),
  FILTER_ID_GR = (1 << 6),
  FILTER_ID_IM = (1 << 7),
  FILTER_ID_LA = (1 << 8),
  FILTER_ID_LS = (1 << 9),
  FILTER_ID_LT = (1 << 10),
  FILTER_ID_MA = (1 << 11),
  FILTER_ID_MB = (1 << 12),
  FILTER_ID_MC = (1 << 13),
  FILTER_ID_ME = (1 << 14),
  FILTER_ID_MSK = (1 << 15),
  FILTER_ID_NT = (1 << 16),
  FILTER_ID_OB = (1 << 17),
  FILTER_ID_PAL = (1 << 18),
  FILTER_ID_PC = (1 << 19),
  FILTER_ID_SCE = (1 << 20),
  FILTER_ID_SPK = (1 << 21),
  FILTER_ID_SO = (1 << 22),
  FILTER_ID_TE = (1 << 23),
  FILTER_ID_TXT = (1 << 24),
  FILTER_ID_VF = (1 << 25),
  FILTER_ID_WO = (1 << 26),
  FILTER_ID_PA = (1 << 27),
  FILTER_ID_CF = (1 << 28),
  FILTER_ID_WS = (1 << 29),
  FILTER_ID_LP = (1u << 31),
};

/* IMPORTANT: this enum matches the order currently use in set_listbasepointers,
 * keep them in sync! */
enum {
  INDEX_ID_LI = 0,
  INDEX_ID_IP,
  INDEX_ID_AC,
  INDEX_ID_KE,
  INDEX_ID_PAL,
  INDEX_ID_GD,
  INDEX_ID_NT,
  INDEX_ID_IM,
  INDEX_ID_TE,
  INDEX_ID_MA,
  INDEX_ID_VF,
  INDEX_ID_AR,
  INDEX_ID_CF,
  INDEX_ID_ME,
  INDEX_ID_CU,
  INDEX_ID_MB,
  INDEX_ID_LT,
  INDEX_ID_LA,
  INDEX_ID_CA,
  INDEX_ID_TXT,
  INDEX_ID_SO,
  INDEX_ID_GR,
  INDEX_ID_PC,
  INDEX_ID_BR,
  INDEX_ID_PA,
  INDEX_ID_SPK,
  INDEX_ID_LP,
  INDEX_ID_WO,
  INDEX_ID_MC,
  INDEX_ID_SCR,
  INDEX_ID_OB,
  INDEX_ID_LS,
  INDEX_ID_SCE,
  INDEX_ID_WS,
  INDEX_ID_WM,
  INDEX_ID_MSK,
  INDEX_ID_NULL,
  INDEX_ID_MAX,
};

#ifdef __cplusplus
}
#endif

#endif

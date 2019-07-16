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
#ifndef __BKE_LIBRARY_H__
#define __BKE_LIBRARY_H__

/** \file
 * \ingroup bke
 */
#ifdef __cplusplus
extern "C" {
#endif

#include "BLI_compiler_attrs.h"

/**
 * Naming: BKE_id_ vs BKE_libblock_
 *
 * WARNING: description below is ideal goal, current status of naming does not yet
 * fully follow it (this is WIP).
 *
 * - BKE_id_ should be used for rather high-level operations, that involve Main database and
 *   relations with other IDs, and can be considered as 'safe'
 *   (as in, in themselves, they leave affected IDs/Main in a consistent status).
 *
 * - BKE_libblock_ should be used for lower level operations,
 *   that perform some parts of BKE_id_ ones, but will generally not ensure caller that affected
 *   data is in a consistent state by their own execution alone.
 *
 * Consequently, external code should not typically use BKE_libblock_ functions,
 * except in some specific cases requiring advanced (and potentially dangerous) handling.
 */

struct BlendThumbnail;
struct GHash;
struct ID;
struct ImBuf;
struct Library;
struct ListBase;
struct Main;
struct PointerRNA;
struct PropertyRNA;
struct bContext;
struct wmWindowManager;

size_t BKE_libblock_get_alloc_info(short type, const char **name);
void *BKE_libblock_alloc_notest(short type) ATTR_WARN_UNUSED_RESULT;
void *BKE_libblock_alloc(struct Main *bmain, short type, const char *name, const int flag)
    ATTR_WARN_UNUSED_RESULT;
void BKE_libblock_init_empty(struct ID *id) ATTR_NONNULL(1);

void *BKE_id_new(struct Main *bmain, const short type, const char *name);
void *BKE_id_new_nomain(const short type, const char *name);

/**
 * New ID creation/copying options.
 */
enum {
  /* *** Generic options (should be handled by all ID types copying, ID creation, etc.). *** */
  /** Create datablock outside of any main database -
   * similar to 'localize' functions of materials etc. */
  LIB_ID_CREATE_NO_MAIN = 1 << 0,
  /** Do not affect user refcount of datablocks used by new one
   * (which also gets zero usercount then).
   * Implies LIB_ID_CREATE_NO_MAIN. */
  LIB_ID_CREATE_NO_USER_REFCOUNT = 1 << 1,
  /** Assume given 'newid' already points to allocated memory for whole datablock
   * (ID + data) - USE WITH CAUTION!
   * Implies LIB_ID_CREATE_NO_MAIN. */
  LIB_ID_CREATE_NO_ALLOCATE = 1 << 2,

  /** Do not tag new ID for update in depsgraph. */
  LIB_ID_CREATE_NO_DEG_TAG = 1 << 8,

  /* *** Specific options to some ID types or usages. *** */
  /* *** May be ignored by unrelated ID copying functions. *** */
  /** Object only, needed by make_local code. */
  LIB_ID_COPY_NO_PROXY_CLEAR = 1 << 16,
  /** Do not copy preview data, when supported. */
  LIB_ID_COPY_NO_PREVIEW = 1 << 17,
  /** Copy runtime data caches. */
  LIB_ID_COPY_CACHES = 1 << 18,
  /** Don't copy id->adt, used by ID datablock localization routines. */
  LIB_ID_COPY_NO_ANIMDATA = 1 << 19,
  /** Mesh: Reference CD data layers instead of doing real copy - USE WITH CAUTION! */
  LIB_ID_COPY_CD_REFERENCE = 1 << 20,

  /* *** XXX Hackish/not-so-nice specific behaviors needed for some corner cases. *** */
  /* *** Ideally we should not have those, but we need them for now... *** */
  /** EXCEPTION! Deep-copy actions used by animdata of copied ID. */
  LIB_ID_COPY_ACTIONS = 1 << 24,
  /** Keep the library pointer when copying datablock outside of bmain. */
  LIB_ID_COPY_KEEP_LIB = 1 << 25,
  /** EXCEPTION! Deep-copy shapekeys used by copied obdata ID. */
  LIB_ID_COPY_SHAPEKEY = 1 << 26,

  /* *** Helper 'defines' gathering most common flag sets. *** */
  /** Shapekeys are not real ID's, more like local data to geometry IDs... */
  LIB_ID_COPY_DEFAULT = LIB_ID_COPY_SHAPEKEY,
  /** Generate a local copy, outside of bmain, to work on (used by COW e.g.). */
  LIB_ID_COPY_LOCALIZE = LIB_ID_CREATE_NO_MAIN | LIB_ID_CREATE_NO_USER_REFCOUNT |
                         LIB_ID_CREATE_NO_DEG_TAG | LIB_ID_COPY_NO_PREVIEW | LIB_ID_COPY_CACHES,
};

void BKE_libblock_copy_ex(struct Main *bmain,
                          const struct ID *id,
                          struct ID **r_newid,
                          const int flag);
void *BKE_libblock_copy(struct Main *bmain, const struct ID *id) ATTR_WARN_UNUSED_RESULT
    ATTR_NONNULL();
/* Special version. sued by datablock localization. */
void *BKE_libblock_copy_for_localize(const struct ID *id);

void BKE_libblock_rename(struct Main *bmain, struct ID *id, const char *name) ATTR_NONNULL();
void BLI_libblock_ensure_unique_name(struct Main *bmain, const char *name) ATTR_NONNULL();

struct ID *BKE_libblock_find_name(struct Main *bmain,
                                  const short type,
                                  const char *name) ATTR_WARN_UNUSED_RESULT ATTR_NONNULL();

/* library_remap.c (keep here since they're general functions) */
/**
 * New freeing logic options.
 */
enum {
  /* *** Generic options (should be handled by all ID types freeing). *** */
  /** Do not try to remove freed ID from given Main (passed Main may be NULL). */
  LIB_ID_FREE_NO_MAIN = 1 << 0,
  /**
   * Do not affect user refcount of datablocks used by freed one.
   * Implies LIB_ID_FREE_NO_MAIN.
   */
  LIB_ID_FREE_NO_USER_REFCOUNT = 1 << 1,
  /**
   * Assume freed ID datablock memory is managed elsewhere, do not free it
   * (still calls relevant ID type's freeing function though) - USE WITH CAUTION!
   * Implies LIB_ID_FREE_NO_MAIN.
   */
  LIB_ID_FREE_NOT_ALLOCATED = 1 << 2,

  /** Do not tag freed ID for update in depsgraph. */
  LIB_ID_FREE_NO_DEG_TAG = 1 << 8,
  /** Do not attempt to remove freed ID from UI data/notifiers/... */
  LIB_ID_FREE_NO_UI_USER = 1 << 9,
};

void BKE_libblock_free_datablock(struct ID *id, const int flag) ATTR_NONNULL();
void BKE_libblock_free_data(struct ID *id, const bool do_id_user) ATTR_NONNULL();

void BKE_id_free_ex(struct Main *bmain, void *idv, int flag, const bool use_flag_from_idtag);
void BKE_id_free(struct Main *bmain, void *idv);

void BKE_id_free_us(struct Main *bmain, void *idv) ATTR_NONNULL();

void BKE_id_delete(struct Main *bmain, void *idv) ATTR_NONNULL();
void BKE_id_multi_tagged_delete(struct Main *bmain) ATTR_NONNULL();

void BKE_libblock_management_main_add(struct Main *bmain, void *idv);
void BKE_libblock_management_main_remove(struct Main *bmain, void *idv);

void BKE_libblock_management_usercounts_set(struct Main *bmain, void *idv);
void BKE_libblock_management_usercounts_clear(struct Main *bmain, void *idv);

void BKE_id_lib_local_paths(struct Main *bmain, struct Library *lib, struct ID *id);
void id_lib_extern(struct ID *id);
void BKE_library_filepath_set(struct Main *bmain, struct Library *lib, const char *filepath);
void id_us_ensure_real(struct ID *id);
void id_us_clear_real(struct ID *id);
void id_us_plus_no_lib(struct ID *id);
void id_us_plus(struct ID *id);
void id_us_min(struct ID *id);
void id_fake_user_set(struct ID *id);
void id_fake_user_clear(struct ID *id);
void BKE_id_clear_newpoin(struct ID *id);

void BKE_id_make_local_generic(struct Main *bmain,
                               struct ID *id,
                               const bool id_in_mainlist,
                               const bool lib_local);
bool id_make_local(struct Main *bmain, struct ID *id, const bool test, const bool force_local);
bool id_single_user(struct bContext *C,
                    struct ID *id,
                    struct PointerRNA *ptr,
                    struct PropertyRNA *prop);
bool BKE_id_copy_is_allowed(const struct ID *id);
bool BKE_id_copy(struct Main *bmain, const struct ID *id, struct ID **newid);
bool BKE_id_copy_ex(struct Main *bmain, const struct ID *id, struct ID **r_newid, const int flag);
void BKE_id_swap(struct Main *bmain, struct ID *id_a, struct ID *id_b);
void id_sort_by_name(struct ListBase *lb, struct ID *id);
void BKE_id_expand_local(struct Main *bmain, struct ID *id);
void BKE_id_copy_ensure_local(struct Main *bmain, const struct ID *old_id, struct ID *new_id);

bool BKE_id_new_name_validate(struct ListBase *lb, struct ID *id, const char *name)
    ATTR_NONNULL(1, 2);
void id_clear_lib_data(struct Main *bmain, struct ID *id);
void id_clear_lib_data_ex(struct Main *bmain, struct ID *id, const bool id_in_mainlist);

/* Affect whole Main database. */
void BKE_main_id_tag_idcode(struct Main *mainvar,
                            const short type,
                            const int tag,
                            const bool value);
void BKE_main_id_tag_listbase(struct ListBase *lb, const int tag, const bool value);
void BKE_main_id_tag_all(struct Main *mainvar, const int tag, const bool value);

void BKE_main_id_flag_listbase(struct ListBase *lb, const int flag, const bool value);
void BKE_main_id_flag_all(struct Main *bmain, const int flag, const bool value);

void BKE_main_id_clear_newpoins(struct Main *bmain);

void BLE_main_id_refcount_recompute(struct Main *bmain, const bool do_linked_only);

void BKE_main_lib_objects_recalc_all(struct Main *bmain);

/* Only for repairing files via versioning, avoid for general use. */
void BKE_main_id_repair_duplicate_names_listbase(struct ListBase *lb);

#define MAX_ID_FULL_NAME (64 + 64 + 3 + 1)         /* 64 is MAX_ID_NAME - 2 */
#define MAX_ID_FULL_NAME_UI (MAX_ID_FULL_NAME + 3) /* Adds 'keycode' two letters at beginning. */
void BKE_id_full_name_get(char name[MAX_ID_FULL_NAME], const struct ID *id);
void BKE_id_full_name_ui_prefix_get(char name[MAX_ID_FULL_NAME_UI], const struct ID *id);

char *BKE_id_to_unique_string_key(const struct ID *id);

void BKE_library_free(struct Library *lib);

void BKE_library_make_local(struct Main *bmain,
                            const struct Library *lib,
                            struct GHash *old_to_new_ids,
                            const bool untagged_only,
                            const bool set_fake);

void BKE_id_tag_set_atomic(struct ID *id, int tag);
void BKE_id_tag_clear_atomic(struct ID *id, int tag);

bool BKE_id_is_in_global_main(struct ID *id);

void BKE_id_ordered_list(struct ListBase *ordered_lb, const struct ListBase *lb);
void BKE_id_reorder(const struct ListBase *lb, struct ID *id, struct ID *relative, bool after);

#define IS_TAGGED(_id) ((_id) && (((ID *)_id)->tag & LIB_TAG_DOIT))

#ifdef __cplusplus
}
#endif

#endif /* __BKE_LIBRARY_H__ */

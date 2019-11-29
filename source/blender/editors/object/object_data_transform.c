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
 * \ingroup edobj
 *
 * Use to transform object origins only.
 *
 * This is a small API to store & apply transformations to object data,
 * where a transformation matrix can be continually applied ontop of the original values
 * so we don't loose precision over time.
 */

#include <stdlib.h>
#include <string.h>

#include "DNA_anim_types.h"
#include "DNA_armature_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meta_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"
#include "DNA_collection_types.h"
#include "DNA_lattice_types.h"

#include "BLI_math.h"
#include "BLI_listbase.h"
#include "BLI_utildefines.h"

#include "BKE_curve.h"
#include "BKE_editmesh.h"
#include "BKE_mball.h"
#include "BKE_mesh.h"
#include "BKE_armature.h"
#include "BKE_lattice.h"

#include "bmesh.h"

#include "DEG_depsgraph.h"

#include "WM_types.h"

#include "ED_object.h"
#include "ED_mesh.h"
#include "ED_armature.h"

#include "MEM_guardedalloc.h"

/* -------------------------------------------------------------------- */
/** \name Internal Transform Get/Apply
 *
 * Some object data types don't have utility functions to access their transformation data.
 * Define these locally.
 *
 * \{ */

/* Armature */

struct ElemData_Armature {
  float tail[3];
  float head[3];
  float roll;
  float arm_tail[3];
  float arm_head[3];
  float arm_roll;
  float rad_tail;
  float rad_head;
  float dist;
  float xwidth;
  float zwidth;
};

static struct ElemData_Armature *armature_coords_and_quats_get_recurse(
    const ListBase *bone_base, struct ElemData_Armature *elem_array)
{
  struct ElemData_Armature *elem = elem_array;
  for (const Bone *bone = bone_base->first; bone; bone = bone->next) {

#define COPY_PTR(member) memcpy(elem->member, bone->member, sizeof(bone->member))
#define COPY_VAL(member) memcpy(&elem->member, &bone->member, sizeof(bone->member))
    COPY_PTR(head);
    COPY_PTR(tail);
    COPY_VAL(roll);
    COPY_PTR(arm_head);
    COPY_PTR(arm_tail);
    COPY_VAL(arm_roll);
    COPY_VAL(rad_tail);
    COPY_VAL(rad_head);
    COPY_VAL(dist);
    COPY_VAL(xwidth);
    COPY_VAL(zwidth);
#undef COPY_PTR
#undef COPY_VAL

    elem = armature_coords_and_quats_get_recurse(&bone->childbase, elem + 1);
  }
  return elem;
}

static void armature_coords_and_quats_get(const bArmature *arm,
                                          struct ElemData_Armature *elem_array)
{
  armature_coords_and_quats_get_recurse(&arm->bonebase, elem_array);
}

static const struct ElemData_Armature *armature_coords_and_quats_apply_with_mat4_recurse(
    ListBase *bone_base, const struct ElemData_Armature *elem_array, const float mat[4][4])
{
  const struct ElemData_Armature *elem = elem_array;
  for (Bone *bone = bone_base->first; bone; bone = bone->next) {

#define COPY_PTR(member) memcpy(bone->member, elem->member, sizeof(bone->member))
#define COPY_VAL(member) memcpy(&bone->member, &elem->member, sizeof(bone->member))
    COPY_PTR(head);
    COPY_PTR(tail);
    COPY_VAL(roll);
    COPY_PTR(arm_head);
    COPY_PTR(arm_tail);
    COPY_VAL(arm_roll);
    COPY_VAL(rad_tail);
    COPY_VAL(rad_head);
    COPY_VAL(dist);
    COPY_VAL(xwidth);
    COPY_VAL(zwidth);
#undef COPY_PTR
#undef COPY_VAL

    elem = armature_coords_and_quats_apply_with_mat4_recurse(&bone->childbase, elem + 1, mat);
  }
  return elem;
}

static void armature_coords_and_quats_apply_with_mat4(bArmature *arm,
                                                      const struct ElemData_Armature *elem_array,
                                                      const float mat[4][4])
{
  armature_coords_and_quats_apply_with_mat4_recurse(&arm->bonebase, elem_array, mat);
  BKE_armature_transform(arm, mat, true);
}

static void armature_coords_and_quats_apply(bArmature *arm,
                                            const struct ElemData_Armature *elem_array)
{
  /* Avoid code duplication by using a unit matrix. */
  float mat[4][4];
  unit_m4(mat);
  armature_coords_and_quats_apply_with_mat4(arm, elem_array, mat);
}

/* Edit Armature */
static void edit_armature_coords_and_quats_get(const bArmature *arm,
                                               struct ElemData_Armature *elem_array)
{
  struct ElemData_Armature *elem = elem_array;
  for (EditBone *ebone = arm->edbo->first; ebone; ebone = ebone->next, elem++) {

#define COPY_PTR(member) memcpy(elem->member, ebone->member, sizeof(ebone->member))
#define COPY_VAL(member) memcpy(&elem->member, &ebone->member, sizeof(ebone->member))
    /* Unused for edit bones: arm_head, arm_tail, arm_roll */
    COPY_PTR(head);
    COPY_PTR(tail);
    COPY_VAL(roll);
    COPY_VAL(rad_tail);
    COPY_VAL(rad_head);
    COPY_VAL(dist);
    COPY_VAL(xwidth);
    COPY_VAL(zwidth);
#undef COPY_PTR
#undef COPY_VAL
  }
}

static void edit_armature_coords_and_quats_apply_with_mat4(
    bArmature *arm, const struct ElemData_Armature *elem_array, const float mat[4][4])
{
  const struct ElemData_Armature *elem = elem_array;
  for (EditBone *ebone = arm->edbo->first; ebone; ebone = ebone->next, elem++) {

#define COPY_PTR(member) memcpy(ebone->member, elem->member, sizeof(ebone->member))
#define COPY_VAL(member) memcpy(&ebone->member, &elem->member, sizeof(ebone->member))
    /* Unused for edit bones: arm_head, arm_tail, arm_roll */
    COPY_PTR(head);
    COPY_PTR(tail);
    COPY_VAL(roll);
    COPY_VAL(rad_tail);
    COPY_VAL(rad_head);
    COPY_VAL(dist);
    COPY_VAL(xwidth);
    COPY_VAL(zwidth);
#undef COPY_PTR
#undef COPY_VAL
  }
  ED_armature_edit_transform(arm, mat, true);
}

static void edit_armature_coords_and_quats_apply(bArmature *arm,
                                                 const struct ElemData_Armature *elem_array)
{
  /* Avoid code duplication by using a unit matrix. */
  float mat[4][4];
  unit_m4(mat);
  edit_armature_coords_and_quats_apply_with_mat4(arm, elem_array, mat);
}

/* MetaBall */

struct ElemData_MetaBall {
  float co[3];
  float quat[4];
  float exp[3];
  float rad;
};

static void metaball_coords_and_quats_get(const MetaBall *mb, struct ElemData_MetaBall *elem_array)
{
  struct ElemData_MetaBall *elem = elem_array;
  for (const MetaElem *ml = mb->elems.first; ml; ml = ml->next, elem++) {
    copy_v3_v3(elem->co, &ml->x);
    copy_qt_qt(elem->quat, ml->quat);
    copy_v3_v3(elem->exp, &ml->expx);
    elem->rad = ml->rad;
  }
}

static void metaball_coords_and_quats_apply_with_mat4(MetaBall *mb,
                                                      const struct ElemData_MetaBall *elem_array,
                                                      const float mat[4][4])
{
  const struct ElemData_MetaBall *elem = elem_array;
  for (MetaElem *ml = mb->elems.first; ml; ml = ml->next, elem++) {
    copy_v3_v3(&ml->x, elem->co);
    copy_qt_qt(ml->quat, elem->quat);
    copy_v3_v3(&ml->expx, elem->exp);
    ml->rad = elem->rad;
  }
  BKE_mball_transform(mb, mat, true);
}

static void metaball_coords_and_quats_apply(MetaBall *mb,
                                            const struct ElemData_MetaBall *elem_array)
{
  /* Avoid code duplication by using a unit matrix. */
  float mat[4][4];
  unit_m4(mat);
  metaball_coords_and_quats_apply_with_mat4(mb, elem_array, mat);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Public Object Data Storage API
 *
 * Used for interactively transforming object data.
 *
 * Store object data transformation in an opaque struct.
 * \{ */

struct XFormObjectData {
  ID *id;
  bool is_edit_mode;
};

struct XFormObjectData_Mesh {
  struct XFormObjectData base;
  float elem_array[0][3];
};

struct XFormObjectData_Lattice {
  struct XFormObjectData base;
  float elem_array[0][3];
};

struct XFormObjectData_Curve {
  struct XFormObjectData base;
  float elem_array[0][3];
};

struct XFormObjectData_Armature {
  struct XFormObjectData base;
  struct ElemData_Armature elem_array[0];
};

struct XFormObjectData_MetaBall {
  struct XFormObjectData base;
  struct ElemData_MetaBall elem_array[0];
};

struct XFormObjectData *ED_object_data_xform_create_ex(ID *id, bool is_edit_mode)
{
  struct XFormObjectData *xod_base = NULL;
  switch (GS(id->name)) {
    case ID_ME: {
      Mesh *me = (Mesh *)id;
      if (is_edit_mode) {
        BMesh *bm = me->edit_mesh->bm;
        const int elem_array_len = bm->totvert;
        struct XFormObjectData_Mesh *xod = MEM_mallocN(
            sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
        BM_mesh_vert_coords_get(bm, xod->elem_array);
        xod_base = &xod->base;
      }
      else {
        const int elem_array_len = me->totvert;
        struct XFormObjectData_Mesh *xod = MEM_mallocN(
            sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
        BKE_mesh_vert_coords_get(me, xod->elem_array);
        xod_base = &xod->base;
      }
      break;
    }
    case ID_LT: {
      Lattice *lt_orig = (Lattice *)id;
      Lattice *lt = is_edit_mode ? lt_orig->editlatt->latt : lt_orig;
      const int elem_array_len = lt->pntsu * lt->pntsv * lt->pntsw;
      struct XFormObjectData_Lattice *xod = MEM_mallocN(
          sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
      BKE_lattice_vert_coords_get(lt, xod->elem_array);
      xod_base = &xod->base;
      break;
    }
    case ID_CU: {
      Curve *cu = (Curve *)id;
      const short ob_type = BKE_curve_type_get(cu);
      if (ob_type == OB_FONT) {
        /* We could support translation. */
        break;
      }

      ListBase *nurbs;
      if (is_edit_mode) {
        EditNurb *editnurb = cu->editnurb;
        nurbs = &editnurb->nurbs;
      }
      else {
        nurbs = &cu->nurb;
      }

      const int elem_array_len = BKE_nurbList_verts_count(nurbs);
      struct XFormObjectData_Curve *xod = MEM_mallocN(
          sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
      BKE_curve_nurbs_vert_coords_get(nurbs, xod->elem_array, elem_array_len);
      xod_base = &xod->base;
      break;
    }
    case ID_AR: {
      bArmature *arm = (bArmature *)id;
      if (is_edit_mode) {
        const int elem_array_len = BLI_listbase_count(arm->edbo);
        struct XFormObjectData_Armature *xod = MEM_mallocN(
            sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
        edit_armature_coords_and_quats_get(arm, xod->elem_array);
        xod_base = &xod->base;
      }
      else {
        const int elem_array_len = BKE_armature_bonelist_count(&arm->bonebase);
        struct XFormObjectData_Armature *xod = MEM_mallocN(
            sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
        armature_coords_and_quats_get(arm, xod->elem_array);
        xod_base = &xod->base;
      }
      break;
    }
    case ID_MB: {
      /* Edit mode and object mode are shared. */
      MetaBall *mb = (MetaBall *)id;
      const int elem_array_len = BLI_listbase_count(&mb->elems);
      struct XFormObjectData_MetaBall *xod = MEM_mallocN(
          sizeof(*xod) + (sizeof(*xod->elem_array) * elem_array_len), __func__);
      metaball_coords_and_quats_get(mb, xod->elem_array);
      xod_base = &xod->base;
      break;
    }
    default: {
      break;
    }
  }
  if (xod_base) {
    xod_base->id = id;
    xod_base->is_edit_mode = is_edit_mode;
  }
  return xod_base;
}

struct XFormObjectData *ED_object_data_xform_create(ID *id)
{
  return ED_object_data_xform_create_ex(id, false);
}

struct XFormObjectData *ED_object_data_xform_create_from_edit_mode(ID *id)
{
  return ED_object_data_xform_create_ex(id, true);
}

void ED_object_data_xform_destroy(struct XFormObjectData *xod)
{
  MEM_freeN(xod);
}

void ED_object_data_xform_by_mat4(struct XFormObjectData *xod_base, const float mat[4][4])
{
  switch (GS(xod_base->id->name)) {
    case ID_ME: {
      Mesh *me = (Mesh *)xod_base->id;
      struct XFormObjectData_Mesh *xod = (struct XFormObjectData_Mesh *)xod_base;
      if (xod_base->is_edit_mode) {
        BMesh *bm = me->edit_mesh->bm;
        BM_mesh_vert_coords_apply_with_mat4(bm, xod->elem_array, mat);
      }
      else {
        BKE_mesh_vert_coords_apply_with_mat4(me, xod->elem_array, mat);
      }
      break;
    }
    case ID_LT: {
      Lattice *lt_orig = (Lattice *)xod_base->id;
      Lattice *lt = xod_base->is_edit_mode ? lt_orig->editlatt->latt : lt_orig;
      struct XFormObjectData_Lattice *xod = (struct XFormObjectData_Lattice *)xod_base;
      BKE_lattice_vert_coords_apply_with_mat4(lt, xod->elem_array, mat);
      break;
    }
    case ID_CU: {
      BLI_assert(xod_base->is_edit_mode == false); /* Not used currently. */
      Curve *cu = (Curve *)xod_base->id;
      struct XFormObjectData_Curve *xod = (struct XFormObjectData_Curve *)xod_base;
      if (xod_base->is_edit_mode) {
        EditNurb *editnurb = cu->editnurb;
        BKE_curve_nurbs_vert_coords_apply_with_mat4(&editnurb->nurbs, xod->elem_array, mat, true);
      }
      else {
        BKE_curve_nurbs_vert_coords_apply_with_mat4(&cu->nurb, xod->elem_array, mat, true);
      }
      break;
    }
    case ID_AR: {
      BLI_assert(xod_base->is_edit_mode == false); /* Not used currently. */
      bArmature *arm = (bArmature *)xod_base->id;
      struct XFormObjectData_Armature *xod = (struct XFormObjectData_Armature *)xod_base;
      if (xod_base->is_edit_mode) {
        edit_armature_coords_and_quats_apply_with_mat4(arm, xod->elem_array, mat);
      }
      else {
        armature_coords_and_quats_apply_with_mat4(arm, xod->elem_array, mat);
      }
      break;
    }
    case ID_MB: {
      /* Metaballs are a special case, edit-mode and object mode data is shared. */
      MetaBall *mb = (MetaBall *)xod_base->id;
      struct XFormObjectData_MetaBall *xod = (struct XFormObjectData_MetaBall *)xod_base;
      metaball_coords_and_quats_apply_with_mat4(mb, xod->elem_array, mat);
      break;
    }
    default: {
      break;
    }
  }
}

void ED_object_data_xform_restore(struct XFormObjectData *xod_base)
{
  switch (GS(xod_base->id->name)) {
    case ID_ME: {
      Mesh *me = (Mesh *)xod_base->id;
      struct XFormObjectData_Mesh *xod = (struct XFormObjectData_Mesh *)xod_base;
      if (xod_base->is_edit_mode) {
        BMesh *bm = me->edit_mesh->bm;
        BM_mesh_vert_coords_apply(bm, xod->elem_array);
      }
      else {
        BKE_mesh_vert_coords_apply(me, xod->elem_array);
      }
      break;
    }
    case ID_LT: {
      Lattice *lt_orig = (Lattice *)xod_base->id;
      Lattice *lt = xod_base->is_edit_mode ? lt_orig->editlatt->latt : lt_orig;
      struct XFormObjectData_Lattice *xod = (struct XFormObjectData_Lattice *)xod_base;
      BKE_lattice_vert_coords_apply(lt, xod->elem_array);
      break;
    }
    case ID_CU: {
      Curve *cu = (Curve *)xod_base->id;
      struct XFormObjectData_Curve *xod = (struct XFormObjectData_Curve *)xod_base;
      if (xod_base->is_edit_mode) {
        EditNurb *editnurb = cu->editnurb;
        BKE_curve_nurbs_vert_coords_apply(&editnurb->nurbs, xod->elem_array, true);
      }
      else {
        BKE_curve_nurbs_vert_coords_apply(&cu->nurb, xod->elem_array, true);
      }
      break;
    }
    case ID_AR: {
      bArmature *arm = (bArmature *)xod_base->id;
      struct XFormObjectData_Armature *xod = (struct XFormObjectData_Armature *)xod_base;
      if (xod_base->is_edit_mode) {
        edit_armature_coords_and_quats_apply(arm, xod->elem_array);
      }
      else {
        armature_coords_and_quats_apply(arm, xod->elem_array);
      }
      break;
    }
    case ID_MB: {
      /* Metaballs are a special case, edit-mode and object mode data is shared. */
      MetaBall *mb = (MetaBall *)xod_base->id;
      struct XFormObjectData_MetaBall *xod = (struct XFormObjectData_MetaBall *)xod_base;
      metaball_coords_and_quats_apply(mb, xod->elem_array);
      break;
    }
    default: {
      break;
    }
  }
}

void ED_object_data_xform_tag_update(struct XFormObjectData *xod_base)
{
  switch (GS(xod_base->id->name)) {
    case ID_ME: {
      Mesh *me = (Mesh *)xod_base->id;
      if (xod_base->is_edit_mode) {
        EDBM_update_generic(me->edit_mesh, true, false);
        EDBM_mesh_normals_update(me->edit_mesh);
      }
      DEG_id_tag_update(&me->id, ID_RECALC_GEOMETRY);
      break;
    }
    case ID_LT: {
      /* Generic update. */
      Lattice *lt = (Lattice *)xod_base->id;
      DEG_id_tag_update(&lt->id, ID_RECALC_GEOMETRY);
      break;
    }
    case ID_CU: {
      /* Generic update. */
      Curve *cu = (Curve *)xod_base->id;
      DEG_id_tag_update(&cu->id, ID_RECALC_GEOMETRY);
      break;
    }
    case ID_AR: {
      /* Generic update. */
      bArmature *arm = (bArmature *)xod_base->id;
      /* XXX, zero is needed, no other flags properly update this. */
      DEG_id_tag_update(&arm->id, 0);
      break;
    }
    case ID_MB: {
      /* Generic update. */
      MetaBall *mb = (MetaBall *)xod_base->id;
      DEG_id_tag_update(&mb->id, ID_RECALC_GEOMETRY | ID_RECALC_COPY_ON_WRITE);
      break;
    }

    default: {
      break;
    }
  }
}

/** \} */

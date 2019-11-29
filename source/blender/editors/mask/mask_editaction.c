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
 * The Original Code is Copyright (C) 2008, Blender Foundation
 * This is a new part of Blender
 */

/** \file
 * \ingroup edmask
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_utildefines.h"

#include "DNA_mask_types.h"
#include "DNA_scene_types.h"

#include "BKE_fcurve.h"
#include "BKE_mask.h"

#include "ED_anim_api.h"
#include "ED_keyframes_edit.h"
#include "ED_mask.h" /* own include */
#include "ED_markers.h"

/* ***************************************** */
/* NOTE ABOUT THIS FILE:
 * This file contains code for editing Mask data in the Action Editor
 * as a 'keyframes', so that a user can adjust the timing of Mask shapekeys.
 * Therefore, this file mostly contains functions for selecting Mask frames (shapekeys).
 */
/* ***************************************** */
/* Generics - Loopers */

/* Loops over the mask-frames for a mask-layer, and applies the given callback */
bool ED_masklayer_frames_looper(MaskLayer *mask_layer,
                                Scene *scene,
                                short (*mask_layer_shape_cb)(MaskLayerShape *, Scene *))
{
  MaskLayerShape *mask_layer_shape;

  /* error checker */
  if (mask_layer == NULL) {
    return false;
  }

  /* do loop */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape->next) {
    /* execute callback */
    if (mask_layer_shape_cb(mask_layer_shape, scene)) {
      return true;
    }
  }

  /* nothing to return */
  return false;
}

/* ****************************************** */
/* Data Conversion Tools */

/* make a listing all the mask-frames in a layer as cfraelems */
void ED_masklayer_make_cfra_list(MaskLayer *mask_layer, ListBase *elems, bool onlysel)
{
  MaskLayerShape *mask_layer_shape;
  CfraElem *ce;

  /* error checking */
  if (ELEM(NULL, mask_layer, elems)) {
    return;
  }

  /* loop through mask-frames, adding */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape->next) {
    if ((onlysel == false) || (mask_layer_shape->flag & MASK_SHAPE_SELECT)) {
      ce = MEM_callocN(sizeof(CfraElem), "CfraElem");

      ce->cfra = (float)mask_layer_shape->frame;
      ce->sel = (mask_layer_shape->flag & MASK_SHAPE_SELECT) ? 1 : 0;

      BLI_addtail(elems, ce);
    }
  }
}

/* ***************************************** */
/* Selection Tools */

/* check if one of the frames in this layer is selected */
bool ED_masklayer_frame_select_check(MaskLayer *mask_layer)
{
  MaskLayerShape *mask_layer_shape;

  /* error checking */
  if (mask_layer == NULL) {
    return 0;
  }

  /* stop at the first one found */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape->next) {
    if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
      return 1;
    }
  }

  /* not found */
  return 0;
}

/* helper function - select mask-frame based on SELECT_* mode */
static void mask_layer_shape_select(MaskLayerShape *mask_layer_shape, short select_mode)
{
  if (mask_layer_shape == NULL) {
    return;
  }

  switch (select_mode) {
    case SELECT_ADD:
      mask_layer_shape->flag |= MASK_SHAPE_SELECT;
      break;
    case SELECT_SUBTRACT:
      mask_layer_shape->flag &= ~MASK_SHAPE_SELECT;
      break;
    case SELECT_INVERT:
      mask_layer_shape->flag ^= MASK_SHAPE_SELECT;
      break;
  }
}

/* set all/none/invert select (like above, but with SELECT_* modes) */
void ED_mask_select_frames(MaskLayer *mask_layer, short select_mode)
{
  MaskLayerShape *mask_layer_shape;

  /* error checking */
  if (mask_layer == NULL) {
    return;
  }

  /* handle according to mode */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape->next) {
    mask_layer_shape_select(mask_layer_shape, select_mode);
  }
}

/* set all/none/invert select */
void ED_masklayer_frame_select_set(MaskLayer *mask_layer, short mode)
{
  /* error checking */
  if (mask_layer == NULL) {
    return;
  }

  /* now call the standard function */
  ED_mask_select_frames(mask_layer, mode);
}

/* select the frame in this layer that occurs on this frame (there should only be one at most) */
void ED_mask_select_frame(MaskLayer *mask_layer, int selx, short select_mode)
{
  MaskLayerShape *mask_layer_shape;

  if (mask_layer == NULL) {
    return;
  }

  mask_layer_shape = BKE_mask_layer_shape_find_frame(mask_layer, selx);

  if (mask_layer_shape) {
    mask_layer_shape_select(mask_layer_shape, select_mode);
  }
}

/* select the frames in this layer that occur within the bounds specified */
void ED_masklayer_frames_select_box(MaskLayer *mask_layer, float min, float max, short select_mode)
{
  MaskLayerShape *mask_layer_shape;

  if (mask_layer == NULL) {
    return;
  }

  /* only select those frames which are in bounds */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape->next) {
    if (IN_RANGE(mask_layer_shape->frame, min, max)) {
      mask_layer_shape_select(mask_layer_shape, select_mode);
    }
  }
}

/* select the frames in this layer that occur within the lasso/circle region specified */
void ED_masklayer_frames_select_region(KeyframeEditData *ked,
                                       MaskLayer *mask_layer,
                                       short tool,
                                       short select_mode)
{
  MaskLayerShape *mask_layer_shape;

  if (mask_layer == NULL) {
    return;
  }

  /* only select frames which are within the region */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape->next) {
    /* construct a dummy point coordinate to do this testing with */
    float pt[2] = {0};

    pt[0] = mask_layer_shape->frame;
    pt[1] = ked->channel_y;

    /* check the necessary regions */
    if (tool == BEZT_OK_CHANNEL_LASSO) {
      /* Lasso */
      if (keyframe_region_lasso_test(ked->data, pt)) {
        mask_layer_shape_select(mask_layer_shape, select_mode);
      }
    }
    else if (tool == BEZT_OK_CHANNEL_CIRCLE) {
      /* Circle */
      if (keyframe_region_circle_test(ked->data, pt)) {
        mask_layer_shape_select(mask_layer_shape, select_mode);
      }
    }
  }
}

/* ***************************************** */
/* Frame Editing Tools */

/* Delete selected frames */
bool ED_masklayer_frames_delete(MaskLayer *mask_layer)
{
  MaskLayerShape *mask_layer_shape, *mask_layer_shape_next;
  bool changed = false;

  /* error checking */
  if (mask_layer == NULL) {
    return false;
  }

  /* check for frames to delete */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = mask_layer_shape_next) {
    mask_layer_shape_next = mask_layer_shape->next;

    if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
      BKE_mask_layer_shape_unlink(mask_layer, mask_layer_shape);
      changed = true;
    }
  }

  return changed;
}

/* Duplicate selected frames from given mask-layer */
void ED_masklayer_frames_duplicate(MaskLayer *mask_layer)
{
  MaskLayerShape *mask_layer_shape, *gpfn;

  /* error checking */
  if (mask_layer == NULL) {
    return;
  }

  /* duplicate selected frames  */
  for (mask_layer_shape = mask_layer->splines_shapes.first; mask_layer_shape;
       mask_layer_shape = gpfn) {
    gpfn = mask_layer_shape->next;

    /* duplicate this frame */
    if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
      MaskLayerShape *mask_shape_dupe;

      /* duplicate frame, and deselect self */
      mask_shape_dupe = BKE_mask_layer_shape_duplicate(mask_layer_shape);
      mask_layer_shape->flag &= ~MASK_SHAPE_SELECT;

      /* XXX - how to handle duplicate frames? */
      BLI_insertlinkafter(&mask_layer->splines_shapes, mask_layer_shape, mask_shape_dupe);
    }
  }
}

/* -------------------------------------- */
/* Snap Tools */

static short snap_mask_layer_nearest(MaskLayerShape *mask_layer_shape, Scene *UNUSED(scene))
{
  if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
    mask_layer_shape->frame = (int)(floor(mask_layer_shape->frame + 0.5));
  }
  return 0;
}

static short snap_mask_layer_nearestsec(MaskLayerShape *mask_layer_shape, Scene *scene)
{
  float secf = (float)FPS;
  if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
    mask_layer_shape->frame = (int)(floorf(mask_layer_shape->frame / secf + 0.5f) * secf);
  }
  return 0;
}

static short snap_mask_layer_cframe(MaskLayerShape *mask_layer_shape, Scene *scene)
{
  if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
    mask_layer_shape->frame = (int)CFRA;
  }
  return 0;
}

static short snap_mask_layer_nearmarker(MaskLayerShape *mask_layer_shape, Scene *scene)
{
  if (mask_layer_shape->flag & MASK_SHAPE_SELECT) {
    mask_layer_shape->frame = (int)ED_markers_find_nearest_marker_time(
        &scene->markers, (float)mask_layer_shape->frame);
  }
  return 0;
}

/* snap selected frames to ... */
void ED_masklayer_snap_frames(MaskLayer *mask_layer, Scene *scene, short mode)
{
  switch (mode) {
    case SNAP_KEYS_NEARFRAME: /* snap to nearest frame */
      ED_masklayer_frames_looper(mask_layer, scene, snap_mask_layer_nearest);
      break;
    case SNAP_KEYS_CURFRAME: /* snap to current frame */
      ED_masklayer_frames_looper(mask_layer, scene, snap_mask_layer_cframe);
      break;
    case SNAP_KEYS_NEARMARKER: /* snap to nearest marker */
      ED_masklayer_frames_looper(mask_layer, scene, snap_mask_layer_nearmarker);
      break;
    case SNAP_KEYS_NEARSEC: /* snap to nearest second */
      ED_masklayer_frames_looper(mask_layer, scene, snap_mask_layer_nearestsec);
      break;
    default: /* just in case */
      break;
  }
}

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
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * Contributor(s): Blender Foundation
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/lattice/editlattice_tools.c
 *  \ingroup edlattice
 */

#include "BLI_math.h"
#include "BLI_utildefines.h"

#include "DNA_curve_types.h"
#include "DNA_lattice_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "BKE_context.h"
#include "BKE_depsgraph.h"
#include "BKE_lattice.h"

#include "ED_screen.h"

#include "WM_api.h"
#include "WM_types.h"

#include "lattice_intern.h"

/** \} */

/* -------------------------------------------------------------------- */
/** \name Make Regular Operator
 * \{ */

static bool make_regular_poll(bContext *C)
{
	Object *ob;

	if (ED_operator_editlattice(C)) return 1;

	ob = CTX_data_active_object(C);
	return (ob && ob->type == OB_LATTICE);
}

static int make_regular_exec(bContext *C, wmOperator *UNUSED(op))
{
	Object *ob = CTX_data_edit_object(C);
	Lattice *lt;

	if (ob) {
		lt = ob->data;
		BKE_lattice_resize(lt->editlatt->latt, lt->pntsu, lt->pntsv, lt->pntsw, NULL);
	}
	else {
		ob = CTX_data_active_object(C);
		lt = ob->data;
		BKE_lattice_resize(lt, lt->pntsu, lt->pntsv, lt->pntsw, NULL);
	}

	DAG_id_tag_update(&ob->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_GEOM | ND_DATA, ob->data);

	return OPERATOR_FINISHED;
}

void LATTICE_OT_make_regular(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Make Regular";
	ot->description = "Set UVW control points a uniform distance apart";
	ot->idname = "LATTICE_OT_make_regular";

	/* api callbacks */
	ot->exec = make_regular_exec;
	ot->poll = make_regular_poll;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Flip Verts Operator
 * \{ */

/* flipping options */
typedef enum eLattice_FlipAxes {
	LATTICE_FLIP_U = 0,
	LATTICE_FLIP_V = 1,
	LATTICE_FLIP_W = 2
} eLattice_FlipAxes;

/**
 * Flip midpoint value so that relative distances between midpoint and neighbor-pair is maintained
 * ! Assumes that uvw <=> xyz (i.e. axis-aligned index-axes with coordinate-axes)
 * - Helper for lattice_flip_exec()
 */
static void lattice_flip_point_value(Lattice *lt, int u, int v, int w, float mid, eLattice_FlipAxes axis)
{
	BPoint *bp;
	float diff;

	/* just the point in the middle (unpaired) */
	bp = &lt->def[BKE_lattice_index_from_uvw(lt, u, v, w)];

	/* flip over axis */
	diff = mid - bp->vec[axis];
	bp->vec[axis] = mid + diff;
}

/**
 * Swap pairs of lattice points along a specified axis
 * - Helper for lattice_flip_exec()
 */
static void lattice_swap_point_pairs(Lattice *lt, int u, int v, int w, float mid, eLattice_FlipAxes axis)
{
	BPoint *bpA, *bpB;

	int numU = lt->pntsu;
	int numV = lt->pntsv;
	int numW = lt->pntsw;

	int u0 = u, u1 = u;
	int v0 = v, v1 = v;
	int w0 = w, w1 = w;

	/* get pair index by just overriding the relevant pair-value
	 * - "-1" else buffer overflow
	 */
	switch (axis) {
		case LATTICE_FLIP_U:
			u1 = numU - u - 1;
			break;
		case LATTICE_FLIP_V:
			v1 = numV - v - 1;
			break;
		case LATTICE_FLIP_W:
			w1 = numW - w - 1;
			break;
	}

	/* get points to operate on */
	bpA = &lt->def[BKE_lattice_index_from_uvw(lt, u0, v0, w0)];
	bpB = &lt->def[BKE_lattice_index_from_uvw(lt, u1, v1, w1)];

	/* Swap all coordinates, so that flipped coordinates belong to
	 * the indices on the correct side of the lattice.
	 *
	 *   Coords:  (-2 4) |0| (3 4)   --> (3 4) |0| (-2 4)
	 *   Indices:  (0,L)     (1,R)   --> (0,L)     (1,R)
	 */
	swap_v3_v3(bpA->vec, bpB->vec);

	/* However, we need to mirror the coordinate values on the axis we're dealing with,
	 * otherwise we'd have effectively only rotated the points around. If we don't do this,
	 * we'd just be reimplementing the naive mirroring algorithm, which causes unwanted deforms
	 * such as flipped normals, etc.
	 *
	 *   Coords:  (3 4) |0| (-2 4)  --\
	 *                                 \-> (-3 4) |0| (2 4)
	 *   Indices: (0,L)     (1,R)   -->     (0,L)     (1,R)
	 */
	lattice_flip_point_value(lt, u0, v0, w0, mid, axis);
	lattice_flip_point_value(lt, u1, v1, w1, mid, axis);
}

static int lattice_flip_exec(bContext *C, wmOperator *op)
{
	Object *obedit = CTX_data_edit_object(C);
	Lattice *lt;

	eLattice_FlipAxes axis = RNA_enum_get(op->ptr, "axis");
	int numU, numV, numW;
	int totP;

	float mid = 0.0f;
	short isOdd = 0;

	/* get lattice - we need the "edit lattice" from the lattice... confusing... */
	lt = (Lattice *)obedit->data;
	lt = lt->editlatt->latt;

	numU = lt->pntsu;
	numV = lt->pntsv;
	numW = lt->pntsw;
	totP = numU * numV * numW;

	/* First Pass: determine midpoint - used for flipping center verts if there are odd number of points on axis */
	switch (axis) {
		case LATTICE_FLIP_U:
			isOdd = numU & 1;
			break;
		case LATTICE_FLIP_V:
			isOdd = numV & 1;
			break;
		case LATTICE_FLIP_W:
			isOdd = numW & 1;
			break;

		default:
			printf("lattice_flip(): Unknown flipping axis (%u)\n", axis);
			return OPERATOR_CANCELLED;
	}

	if (isOdd) {
		BPoint *bp;
		float avgInv = 1.0f / (float)totP;
		int i;

		/* midpoint calculation - assuming that u/v/w are axis-aligned */
		for (i = 0, bp = lt->def; i < totP; i++, bp++) {
			mid += bp->vec[axis] * avgInv;
		}
	}

	/* Second Pass: swap pairs of vertices per axis, assuming they are all sorted */
	switch (axis) {
		case LATTICE_FLIP_U:
		{
			int u, v, w;

			/* v/w strips - front to back, top to bottom */
			for (w = 0; w < numW; w++) {
				for (v = 0; v < numV; v++) {
					/* swap coordinates of pairs of vertices on u */
					for (u = 0; u < (numU / 2); u++) {
						lattice_swap_point_pairs(lt, u, v, w, mid, axis);
					}

					/* flip u-coordinate of midpoint (i.e. unpaired point on u) */
					if (isOdd) {
						u = (numU / 2);
						lattice_flip_point_value(lt, u, v, w, mid, axis);
					}
				}
			}
			break;
		}
		case LATTICE_FLIP_V:
		{
			int u, v, w;

			/* u/w strips - front to back, left to right */
			for (w = 0; w < numW; w++) {
				for (u = 0; u < numU; u++) {
					/* swap coordinates of pairs of vertices on v */
					for (v = 0; v < (numV / 2); v++) {
						lattice_swap_point_pairs(lt, u, v, w, mid, axis);
					}

					/* flip v-coordinate of midpoint (i.e. unpaired point on v) */
					if (isOdd) {
						v = (numV / 2);
						lattice_flip_point_value(lt, u, v, w, mid, axis);
					}
				}
			}
			break;
		}
		case LATTICE_FLIP_W:
		{
			int u, v, w;

			for (v = 0; v < numV; v++) {
				for (u = 0; u < numU; u++) {
					/* swap coordinates of pairs of vertices on w */
					for (w = 0; w < (numW / 2); w++) {
						lattice_swap_point_pairs(lt, u, v, w, mid, axis);
					}

					/* flip w-coordinate of midpoint (i.e. unpaired point on w) */
					if (isOdd) {
						w = (numW / 2);
						lattice_flip_point_value(lt, u, v, w, mid, axis);
					}
				}
			}
			break;
		}
		default: /* shouldn't happen, but just in case */
			break;
	}

	/* updates */
	DAG_id_tag_update(&obedit->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_GEOM | ND_DATA, obedit->data);

	return OPERATOR_FINISHED;
}

void LATTICE_OT_flip(wmOperatorType *ot)
{
	static const EnumPropertyItem flip_items[] = {
		{LATTICE_FLIP_U, "U", 0, "U (X) Axis", ""},
		{LATTICE_FLIP_V, "V", 0, "V (Y) Axis", ""},
		{LATTICE_FLIP_W, "W", 0, "W (Z) Axis", ""},
		{0, NULL, 0, NULL, NULL}};

	/* identifiers */
	ot->name = "Flip (Distortion Free)";
	ot->description = "Mirror all control points without inverting the lattice deform";
	ot->idname = "LATTICE_OT_flip";

	/* api callbacks */
	ot->poll = ED_operator_editlattice;
	ot->invoke = WM_menu_invoke;
	ot->exec = lattice_flip_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	/* properties */
	ot->prop = RNA_def_enum(ot->srna, "axis", flip_items, LATTICE_FLIP_U, "Flip Axis", "Coordinates along this axis get flipped");
}

/** \} */

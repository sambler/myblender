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
 * Contributor(s): Geoffrey Bantle.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/bmesh/intern/bmesh_mesh.c
 *  \ingroup bmesh
 *
 * BM mesh level functions.
 */

#include "MEM_guardedalloc.h"

#include "DNA_listBase.h"
#include "DNA_object_types.h"

#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_utildefines.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_editmesh.h"
#include "BKE_multires.h"

#include "intern/bmesh_private.h"

/* used as an extern, defined in bmesh.h */
const BMAllocTemplate bm_mesh_allocsize_default = {512, 1024, 2048, 512};
const BMAllocTemplate bm_mesh_chunksize_default = {512, 1024, 2048, 512};

static void bm_mempool_init(BMesh *bm, const BMAllocTemplate *allocsize)
{
	bm->vpool = BLI_mempool_create(sizeof(BMVert), allocsize->totvert,
	                               bm_mesh_chunksize_default.totvert, BLI_MEMPOOL_ALLOW_ITER);
	bm->epool = BLI_mempool_create(sizeof(BMEdge), allocsize->totedge,
	                               bm_mesh_chunksize_default.totedge, BLI_MEMPOOL_ALLOW_ITER);
	bm->lpool = BLI_mempool_create(sizeof(BMLoop), allocsize->totloop,
	                               bm_mesh_chunksize_default.totloop, BLI_MEMPOOL_NOP);
	bm->fpool = BLI_mempool_create(sizeof(BMFace), allocsize->totface,
	                               bm_mesh_chunksize_default.totface, BLI_MEMPOOL_ALLOW_ITER);

#ifdef USE_BMESH_HOLES
	bm->looplistpool = BLI_mempool_create(sizeof(BMLoopList), 512, 512, BLI_MEMPOOL_NOP);
#endif
}

void BM_mesh_elem_toolflags_ensure(BMesh *bm)
{
	if (bm->vtoolflagpool && bm->etoolflagpool && bm->ftoolflagpool) {
		return;
	}

	bm->vtoolflagpool = BLI_mempool_create(sizeof(BMFlagLayer), bm->totvert, 512, BLI_MEMPOOL_NOP);
	bm->etoolflagpool = BLI_mempool_create(sizeof(BMFlagLayer), bm->totedge, 512, BLI_MEMPOOL_NOP);
	bm->ftoolflagpool = BLI_mempool_create(sizeof(BMFlagLayer), bm->totface, 512, BLI_MEMPOOL_NOP);

#pragma omp parallel sections if (bm->totvert + bm->totedge + bm->totface >= BM_OMP_LIMIT)
	{
#pragma omp section
		{
			BLI_mempool *toolflagpool = bm->vtoolflagpool;
			BMIter iter;
			BMElemF *ele;
			BM_ITER_MESH (ele, &iter, bm, BM_VERTS_OF_MESH) {
				ele->oflags = BLI_mempool_calloc(toolflagpool);
			}
		}
#pragma omp section
		{
			BLI_mempool *toolflagpool = bm->etoolflagpool;
			BMIter iter;
			BMElemF *ele;
			BM_ITER_MESH (ele, &iter, bm, BM_EDGES_OF_MESH) {
				ele->oflags = BLI_mempool_calloc(toolflagpool);
			}
		}
#pragma omp section
		{
			BLI_mempool *toolflagpool = bm->ftoolflagpool;
			BMIter iter;
			BMElemF *ele;
			BM_ITER_MESH (ele, &iter, bm, BM_FACES_OF_MESH) {
				ele->oflags = BLI_mempool_calloc(toolflagpool);
			}
		}
	}


	bm->totflags = 1;
}

void BM_mesh_elem_toolflags_clear(BMesh *bm)
{
	if (bm->vtoolflagpool) {
		BLI_mempool_destroy(bm->vtoolflagpool);
		bm->vtoolflagpool = NULL;
	}
	if (bm->etoolflagpool) {
		BLI_mempool_destroy(bm->etoolflagpool);
		bm->etoolflagpool = NULL;
	}
	if (bm->ftoolflagpool) {
		BLI_mempool_destroy(bm->ftoolflagpool);
		bm->ftoolflagpool = NULL;
	}
}

/**
 * \brief BMesh Make Mesh
 *
 * Allocates a new BMesh structure.
 *
 * \return The New bmesh
 *
 * \note ob is needed by multires
 */
BMesh *BM_mesh_create(const BMAllocTemplate *allocsize)
{
	/* allocate the structure */
	BMesh *bm = MEM_callocN(sizeof(BMesh), __func__);
	
	/* allocate the memory pools for the mesh elements */
	bm_mempool_init(bm, allocsize);

	/* allocate one flag pool that we don't get rid of. */
	bm->stackdepth = 1;
	bm->totflags = 0;

	CustomData_reset(&bm->vdata);
	CustomData_reset(&bm->edata);
	CustomData_reset(&bm->ldata);
	CustomData_reset(&bm->pdata);

	return bm;
}

/**
 * \brief BMesh Free Mesh Data
 *
 *	Frees a BMesh structure.
 *
 * \note frees mesh, but not actual BMesh struct
 */
void BM_mesh_data_free(BMesh *bm)
{
	BMVert *v;
	BMEdge *e;
	BMLoop *l;
	BMFace *f;

	BMIter iter;
	BMIter itersub;

	const bool is_ldata_free = CustomData_bmesh_has_free(&bm->ldata);
	const bool is_pdata_free = CustomData_bmesh_has_free(&bm->pdata);

	/* Check if we have to call free, if not we can avoid a lot of looping */
	if (CustomData_bmesh_has_free(&(bm->vdata))) {
		BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
			CustomData_bmesh_free_block(&(bm->vdata), &(v->head.data));
		}
	}
	if (CustomData_bmesh_has_free(&(bm->edata))) {
		BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
			CustomData_bmesh_free_block(&(bm->edata), &(e->head.data));
		}
	}

	if (is_ldata_free || is_pdata_free) {
		BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
			if (is_pdata_free)
				CustomData_bmesh_free_block(&(bm->pdata), &(f->head.data));
			if (is_ldata_free) {
				BM_ITER_ELEM (l, &itersub, f, BM_LOOPS_OF_FACE) {
					CustomData_bmesh_free_block(&(bm->ldata), &(l->head.data));
				}
			}
		}
	}

	/* Free custom data pools, This should probably go in CustomData_free? */
	if (bm->vdata.totlayer) BLI_mempool_destroy(bm->vdata.pool);
	if (bm->edata.totlayer) BLI_mempool_destroy(bm->edata.pool);
	if (bm->ldata.totlayer) BLI_mempool_destroy(bm->ldata.pool);
	if (bm->pdata.totlayer) BLI_mempool_destroy(bm->pdata.pool);

	/* free custom data */
	CustomData_free(&bm->vdata, 0);
	CustomData_free(&bm->edata, 0);
	CustomData_free(&bm->ldata, 0);
	CustomData_free(&bm->pdata, 0);

	/* destroy element pools */
	BLI_mempool_destroy(bm->vpool);
	BLI_mempool_destroy(bm->epool);
	BLI_mempool_destroy(bm->lpool);
	BLI_mempool_destroy(bm->fpool);

	if (bm->vtable) MEM_freeN(bm->vtable);
	if (bm->etable) MEM_freeN(bm->etable);
	if (bm->ftable) MEM_freeN(bm->ftable);

	/* destroy flag pool */
	BM_mesh_elem_toolflags_clear(bm);

#ifdef USE_BMESH_HOLES
	BLI_mempool_destroy(bm->looplistpool);
#endif

	BLI_freelistN(&bm->selected);

	BMO_error_clear(bm);
}

/**
 * \brief BMesh Clear Mesh
 *
 * Clear all data in bm
 */
void BM_mesh_clear(BMesh *bm)
{
	/* free old mesh */
	BM_mesh_data_free(bm);
	memset(bm, 0, sizeof(BMesh));

	/* allocate the memory pools for the mesh elements */
	bm_mempool_init(bm, &bm_mesh_allocsize_default);

	bm->stackdepth = 1;
	bm->totflags = 0;

	CustomData_reset(&bm->vdata);
	CustomData_reset(&bm->edata);
	CustomData_reset(&bm->ldata);
	CustomData_reset(&bm->pdata);
}

/**
 * \brief BMesh Free Mesh
 *
 *	Frees a BMesh data and its structure.
 */
void BM_mesh_free(BMesh *bm)
{
	BM_mesh_data_free(bm);

	if (bm->py_handle) {
		/* keep this out of 'BM_mesh_data_free' because we want python
		 * to be able to clear the mesh and maintain access. */
		bpy_bm_generic_invalidate(bm->py_handle);
		bm->py_handle = NULL;
	}

	MEM_freeN(bm);
}

/**
 * Helpers for #BM_mesh_normals_update and #BM_verts_calc_normal_vcos
 */
static void bm_mesh_edges_calc_vectors(BMesh *bm, float (*edgevec)[3], const float (*vcos)[3])
{
	BMIter eiter;
	BMEdge *e;
	int index;

	if (vcos) {
		BM_mesh_elem_index_ensure(bm, BM_VERT);
	}

	BM_ITER_MESH_INDEX (e, &eiter, bm, BM_EDGES_OF_MESH, index) {
		BM_elem_index_set(e, index); /* set_inline */

		if (e->l) {
			const float *v1_co = vcos ? vcos[BM_elem_index_get(e->v1)] : e->v1->co;
			const float *v2_co = vcos ? vcos[BM_elem_index_get(e->v2)] : e->v2->co;
			sub_v3_v3v3(edgevec[index], v2_co, v1_co);
			normalize_v3(edgevec[index]);
		}
		else {
			/* the edge vector will not be needed when the edge has no radial */
		}
	}
	bm->elem_index_dirty &= ~BM_EDGE;
}

static void bm_mesh_verts_calc_normals(BMesh *bm, const float (*edgevec)[3], const float (*fnos)[3],
                                       const float (*vcos)[3], float (*vnos)[3])
{
	BM_mesh_elem_index_ensure(bm, (vnos) ? (BM_EDGE | BM_VERT) : BM_EDGE);

	/* add weighted face normals to vertices */
	{
		BMIter fiter;
		BMFace *f;
		int i;

		BM_ITER_MESH_INDEX (f, &fiter, bm, BM_FACES_OF_MESH, i) {
			BMLoop *l_first, *l_iter;
			const float *f_no = fnos ? fnos[i] : f->no;

			l_iter = l_first = BM_FACE_FIRST_LOOP(f);
			do {
				const float *e1diff, *e2diff;
				float dotprod;
				float fac;
				float *v_no = vnos ? vnos[BM_elem_index_get(l_iter->v)] : l_iter->v->no;

				/* calculate the dot product of the two edges that
				 * meet at the loop's vertex */
				e1diff = edgevec[BM_elem_index_get(l_iter->prev->e)];
				e2diff = edgevec[BM_elem_index_get(l_iter->e)];
				dotprod = dot_v3v3(e1diff, e2diff);

				/* edge vectors are calculated from e->v1 to e->v2, so
				 * adjust the dot product if one but not both loops
				 * actually runs from from e->v2 to e->v1 */
				if ((l_iter->prev->e->v1 == l_iter->prev->v) ^ (l_iter->e->v1 == l_iter->v)) {
					dotprod = -dotprod;
				}

				fac = saacos(-dotprod);

				/* accumulate weighted face normal into the vertex's normal */
				madd_v3_v3fl(v_no, f_no, fac);
			} while ((l_iter = l_iter->next) != l_first);
		}
	}


	/* normalize the accumulated vertex normals */
	{
		BMIter viter;
		BMVert *v;
		int i;

		BM_ITER_MESH_INDEX (v, &viter, bm, BM_VERTS_OF_MESH, i) {
			float *v_no = vnos ? vnos[i] : v->no;
			if (UNLIKELY(normalize_v3(v_no) == 0.0f)) {
				const float *v_co = vcos ? vcos[i] : v->co;
				normalize_v3_v3(v_no, v_co);
			}
		}
	}
}

/**
 * \brief BMesh Compute Normals
 *
 * Updates the normals of a mesh.
 */
void BM_mesh_normals_update(BMesh *bm)
{
	float (*edgevec)[3] = MEM_mallocN(sizeof(*edgevec) * bm->totedge, __func__);

#pragma omp parallel sections if (bm->totvert + bm->totedge + bm->totface >= BM_OMP_LIMIT)
	{
#pragma omp section
		{
			/* calculate all face normals */
			BMIter fiter;
			BMFace *f;
			int i;

			BM_ITER_MESH_INDEX (f, &fiter, bm, BM_FACES_OF_MESH, i) {
				BM_elem_index_set(f, i); /* set_inline */
				BM_face_normal_update(f);
			}
			bm->elem_index_dirty &= ~BM_FACE;
		}
#pragma omp section
		{
			/* Zero out vertex normals */
			BMIter viter;
			BMVert *v;
			int i;

			BM_ITER_MESH_INDEX (v, &viter, bm, BM_VERTS_OF_MESH, i) {
				BM_elem_index_set(v, i); /* set_inline */
				zero_v3(v->no);
			}
			bm->elem_index_dirty &= ~BM_VERT;
		}
#pragma omp section
		{
			/* Compute normalized direction vectors for each edge.
			 * Directions will be used for calculating the weights of the face normals on the vertex normals.
			 */
			bm_mesh_edges_calc_vectors(bm, edgevec, NULL);
		}
	}
	/* end omp */

	/* Add weighted face normals to vertices, and normalize vert normals. */
	bm_mesh_verts_calc_normals(bm, (const float(*)[3])edgevec, NULL, NULL, NULL);
	MEM_freeN(edgevec);
}

/**
 * \brief BMesh Compute Normals from/to external data.
 *
 * Computes the vertex normals of a mesh into vnos, using given vertex coordinates (vcos) and polygon normals (fnos).
 */
void BM_verts_calc_normal_vcos(BMesh *bm, const float (*fnos)[3], const float (*vcos)[3], float (*vnos)[3])
{
	float (*edgevec)[3] = MEM_mallocN(sizeof(*edgevec) * bm->totedge, __func__);

	/* Compute normalized direction vectors for each edge.
	 * Directions will be used for calculating the weights of the face normals on the vertex normals.
	 */
	bm_mesh_edges_calc_vectors(bm, edgevec, vcos);

	/* Add weighted face normals to vertices, and normalize vert normals. */
	bm_mesh_verts_calc_normals(bm, (const float(*)[3])edgevec, fnos, vcos, vnos);
	MEM_freeN(edgevec);
}

static void UNUSED_FUNCTION(bm_mdisps_space_set)(Object *ob, BMesh *bm, int from, int to)
{
	/* switch multires data out of tangent space */
	if (CustomData_has_layer(&bm->ldata, CD_MDISPS)) {
		BMEditMesh *em = BKE_editmesh_create(bm, false);
		DerivedMesh *dm = CDDM_from_editbmesh(em, true, false);
		MDisps *mdisps;
		BMFace *f;
		BMIter iter;
		// int i = 0; // UNUSED
		
		multires_set_space(dm, ob, from, to);
		
		mdisps = CustomData_get_layer(&dm->loopData, CD_MDISPS);
		
		BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
			BMLoop *l;
			BMIter liter;
			BM_ITER_ELEM (l, &liter, f, BM_LOOPS_OF_FACE) {
				MDisps *lmd = CustomData_bmesh_get(&bm->ldata, l->head.data, CD_MDISPS);
				
				if (!lmd->disps) {
					printf("%s: warning - 'lmd->disps' == NULL\n", __func__);
				}
				
				if (lmd->disps && lmd->totdisp == mdisps->totdisp) {
					memcpy(lmd->disps, mdisps->disps, sizeof(float) * 3 * lmd->totdisp);
				}
				else if (mdisps->disps) {
					if (lmd->disps)
						MEM_freeN(lmd->disps);
					
					lmd->disps = MEM_dupallocN(mdisps->disps);
					lmd->totdisp = mdisps->totdisp;
					lmd->level = mdisps->level;
				}
				
				mdisps++;
				// i += 1;
			}
		}
		
		dm->needsFree = 1;
		dm->release(dm);
		
		/* setting this to NULL prevents BKE_editmesh_free from freeing it */
		em->bm = NULL;
		BKE_editmesh_free(em);
		MEM_freeN(em);
	}
}

/**
 * \brief BMesh Begin Edit
 *
 * Functions for setting up a mesh for editing and cleaning up after
 * the editing operations are done. These are called by the tools/operator
 * API for each time a tool is executed.
 */
void bmesh_edit_begin(BMesh *UNUSED(bm), BMOpTypeFlag UNUSED(type_flag))
{
	/* Most operators seem to be using BMO_OPTYPE_FLAG_UNTAN_MULTIRES to change the MDisps to
	 * absolute space during mesh edits. With this enabled, changes to the topology
	 * (loop cuts, edge subdivides, etc) are not reflected in the higher levels of
	 * the mesh at all, which doesn't seem right. Turning off completely for now,
	 * until this is shown to be better for certain types of mesh edits. */
#ifdef BMOP_UNTAN_MULTIRES_ENABLED
	/* switch multires data out of tangent space */
	if ((type_flag & BMO_OPTYPE_FLAG_UNTAN_MULTIRES) && CustomData_has_layer(&bm->ldata, CD_MDISPS)) {
		bmesh_mdisps_space_set(bm, MULTIRES_SPACE_TANGENT, MULTIRES_SPACE_ABSOLUTE);

		/* ensure correct normals, if possible */
		bmesh_rationalize_normals(bm, 0);
		BM_mesh_normals_update(bm);
	}
#endif
}

/**
 * \brief BMesh End Edit
 */
void bmesh_edit_end(BMesh *bm, BMOpTypeFlag type_flag)
{
	/* BMO_OPTYPE_FLAG_UNTAN_MULTIRES disabled for now, see comment above in bmesh_edit_begin. */
#ifdef BMOP_UNTAN_MULTIRES_ENABLED
	/* switch multires data into tangent space */
	if ((flag & BMO_OPTYPE_FLAG_UNTAN_MULTIRES) && CustomData_has_layer(&bm->ldata, CD_MDISPS)) {
		/* set normals to their previous winding */
		bmesh_rationalize_normals(bm, 1);
		bmesh_mdisps_space_set(bm, MULTIRES_SPACE_ABSOLUTE, MULTIRES_SPACE_TANGENT);
	}
	else if (flag & BMO_OP_FLAG_RATIONALIZE_NORMALS) {
		bmesh_rationalize_normals(bm, 1);
	}
#endif

	/* compute normals, clear temp flags and flush selections */
	if (type_flag & BMO_OPTYPE_FLAG_NORMALS_CALC) {
		BM_mesh_normals_update(bm);
	}

	if (type_flag & BMO_OPTYPE_FLAG_SELECT_FLUSH) {
		BM_mesh_select_mode_flush(bm);
	}
}

void BM_mesh_elem_index_ensure(BMesh *bm, const char hflag)
{
#ifdef DEBUG
	BM_ELEM_INDEX_VALIDATE(bm, "Should Never Fail!", __func__);
#endif

#pragma omp parallel sections if (bm->totvert + bm->totedge + bm->totface >= BM_OMP_LIMIT)
	{
#pragma omp section
		{
			if (hflag & BM_VERT) {
				if (bm->elem_index_dirty & BM_VERT) {
					BMIter iter;
					BMElem *ele;

					int index;
					BM_ITER_MESH_INDEX (ele, &iter, bm, BM_VERTS_OF_MESH, index) {
						BM_elem_index_set(ele, index); /* set_ok */
					}
					BLI_assert(index == bm->totvert);
				}
				else {
					// printf("%s: skipping vert index calc!\n", __func__);
				}
			}
		}

#pragma omp section
		{
			if (hflag & BM_EDGE) {
				if (bm->elem_index_dirty & BM_EDGE) {
					BMIter iter;
					BMElem *ele;

					int index;
					BM_ITER_MESH_INDEX (ele, &iter, bm, BM_EDGES_OF_MESH, index) {
						BM_elem_index_set(ele, index); /* set_ok */
					}
					BLI_assert(index == bm->totedge);
				}
				else {
					// printf("%s: skipping edge index calc!\n", __func__);
				}
			}
		}

#pragma omp section
		{
			if (hflag & BM_FACE) {
				if (bm->elem_index_dirty & BM_FACE) {
					BMIter iter;
					BMElem *ele;

					int index;
					BM_ITER_MESH_INDEX (ele, &iter, bm, BM_FACES_OF_MESH, index) {
						BM_elem_index_set(ele, index); /* set_ok */
					}
					BLI_assert(index == bm->totface);
				}
				else {
					// printf("%s: skipping face index calc!\n", __func__);
				}
			}
		}
	}

	bm->elem_index_dirty &= ~hflag;
}


/**
 * Array checking/setting macros
 *
 * Currently vert/edge/loop/face index data is being abused, in a few areas of the code.
 *
 * To avoid correcting them afterwards, set 'bm->elem_index_dirty' however its possible
 * this flag is set incorrectly which could crash blender.
 *
 * These functions ensure its correct and are called more often in debug mode.
 */

void BM_mesh_elem_index_validate(BMesh *bm, const char *location, const char *func,
                                 const char *msg_a, const char *msg_b)
{
	const char iter_types[3] = {BM_VERTS_OF_MESH,
	                            BM_EDGES_OF_MESH,
	                            BM_FACES_OF_MESH};

	const char flag_types[3] = {BM_VERT, BM_EDGE, BM_FACE};
	const char *type_names[3] = {"vert", "edge", "face"};

	BMIter iter;
	BMElem *ele;
	int i;
	bool is_any_error = 0;

	for (i = 0; i < 3; i++) {
		const bool is_dirty = (flag_types[i] & bm->elem_index_dirty) != 0;
		int index = 0;
		bool is_error = false;
		int err_val = 0;
		int err_idx = 0;

		BM_ITER_MESH (ele, &iter, bm, iter_types[i]) {
			if (!is_dirty) {
				if (BM_elem_index_get(ele) != index) {
					err_val = BM_elem_index_get(ele);
					err_idx = index;
					is_error = true;
				}
			}

			BM_elem_index_set(ele, index); /* set_ok */
			index++;
		}

		if ((is_error == true) && (is_dirty == false)) {
			is_any_error = true;
			fprintf(stderr,
			        "Invalid Index: at %s, %s, %s[%d] invalid index %d, '%s', '%s'\n",
			        location, func, type_names[i], err_idx, err_val, msg_a, msg_b);
		}
		else if ((is_error == false) && (is_dirty == true)) {

#if 0       /* mostly annoying */

			/* dirty may have been incorrectly set */
			fprintf(stderr,
			        "Invalid Dirty: at %s, %s (%s), dirty flag was set but all index values are correct, '%s', '%s'\n",
			        location, func, type_names[i], msg_a, msg_b);
#endif
		}
	}

#if 0 /* mostly annoying, even in debug mode */
#ifdef DEBUG
	if (is_any_error == 0) {
		fprintf(stderr,
		        "Valid Index Success: at %s, %s, '%s', '%s'\n",
		        location, func, msg_a, msg_b);
	}
#endif
#endif
	(void) is_any_error; /* shut up the compiler */

}

/* debug check only - no need to optimize */
#ifndef NDEBUG
bool BM_mesh_elem_table_check(BMesh *bm)
{
	BMIter iter;
	BMElem *ele;
	int i;

	if (bm->vtable && ((bm->elem_table_dirty & BM_VERT) == 0)) {
		BM_ITER_MESH_INDEX (ele, &iter, bm, BM_VERTS_OF_MESH, i) {
			if (ele != (BMElem *)bm->vtable[i]) {
				return false;
			}
		}
	}

	if (bm->etable && ((bm->elem_table_dirty & BM_EDGE) == 0)) {
		BM_ITER_MESH_INDEX (ele, &iter, bm, BM_EDGES_OF_MESH, i) {
			if (ele != (BMElem *)bm->etable[i]) {
				return false;
			}
		}
	}

	if (bm->ftable && ((bm->elem_table_dirty & BM_FACE) == 0)) {
		BM_ITER_MESH_INDEX (ele, &iter, bm, BM_FACES_OF_MESH, i) {
			if (ele != (BMElem *)bm->ftable[i]) {
				return false;
			}
		}
	}

	return true;
}
#endif



void BM_mesh_elem_table_ensure(BMesh *bm, const char htype)
{
	/* assume if the array is non-null then its valid and no need to recalc */
	const char htype_needed = (((bm->vtable && ((bm->elem_table_dirty & BM_VERT) == 0)) ? 0 : BM_VERT) |
	                           ((bm->etable && ((bm->elem_table_dirty & BM_EDGE) == 0)) ? 0 : BM_EDGE) |
	                           ((bm->ftable && ((bm->elem_table_dirty & BM_FACE) == 0)) ? 0 : BM_FACE)) & htype;

	BLI_assert((htype & ~BM_ALL_NOLOOP) == 0);

	/* in debug mode double check we didn't need to recalculate */
	BLI_assert(BM_mesh_elem_table_check(bm) == true);

	if (htype_needed & BM_VERT) {
		if (bm->vtable && bm->totvert <= bm->vtable_tot && bm->totvert * 2 >= bm->vtable_tot) {
			/* pass (re-use the array) */
		}
		else {
			if (bm->vtable)
				MEM_freeN(bm->vtable);
			bm->vtable = MEM_mallocN(sizeof(void **) * bm->totvert, "bm->vtable");
			bm->vtable_tot = bm->totvert;
		}
	}
	if (htype_needed & BM_EDGE) {
		if (bm->etable && bm->totedge <= bm->etable_tot && bm->totedge * 2 >= bm->etable_tot) {
			/* pass (re-use the array) */
		}
		else {
			if (bm->etable)
				MEM_freeN(bm->etable);
			bm->etable = MEM_mallocN(sizeof(void **) * bm->totedge, "bm->etable");
			bm->etable_tot = bm->totedge;
		}
	}
	if (htype_needed & BM_FACE) {
		if (bm->ftable && bm->totface <= bm->ftable_tot && bm->totface * 2 >= bm->ftable_tot) {
			/* pass (re-use the array) */
		}
		else {
			if (bm->ftable)
				MEM_freeN(bm->ftable);
			bm->ftable = MEM_mallocN(sizeof(void **) * bm->totface, "bm->ftable");
			bm->ftable_tot = bm->totface;
		}
	}

#pragma omp parallel sections if (bm->totvert + bm->totedge + bm->totface >= BM_OMP_LIMIT)
	{
#pragma omp section
		{
			if (htype_needed & BM_VERT) {
				BM_iter_as_array(bm, BM_VERTS_OF_MESH, NULL, (void **)bm->vtable, bm->totvert);
			}
		}
#pragma omp section
		{
			if (htype_needed & BM_EDGE) {
				BM_iter_as_array(bm, BM_EDGES_OF_MESH, NULL, (void **)bm->etable, bm->totedge);
			}
		}
#pragma omp section
		{
			if (htype_needed & BM_FACE) {
				BM_iter_as_array(bm, BM_FACES_OF_MESH, NULL, (void **)bm->ftable, bm->totface);
			}
		}
	}

	/* Only clear dirty flags when all the pointers and data are actually valid.
	 * This prevents possible threading issues when dirty flag check failed but
	 * data wasn't ready still.
	 */
	if (htype_needed & BM_VERT) {
		bm->elem_table_dirty &= ~BM_VERT;
	}
	if (htype_needed & BM_EDGE) {
		bm->elem_table_dirty &= ~BM_EDGE;
	}
	if (htype_needed & BM_FACE) {
		bm->elem_table_dirty &= ~BM_FACE;
	}
}

/* use BM_mesh_elem_table_ensure where possible to avoid full rebuild */
void BM_mesh_elem_table_init(BMesh *bm, const char htype)
{
	BLI_assert((htype & ~BM_ALL_NOLOOP) == 0);

	/* force recalc */
	BM_mesh_elem_table_free(bm, BM_ALL_NOLOOP);
	BM_mesh_elem_table_ensure(bm, htype);
}

void BM_mesh_elem_table_free(BMesh *bm, const char htype)
{
	if (htype & BM_VERT) {
		MEM_SAFE_FREE(bm->vtable);
	}

	if (htype & BM_EDGE) {
		MEM_SAFE_FREE(bm->etable);
	}

	if (htype & BM_FACE) {
		MEM_SAFE_FREE(bm->ftable);
	}
}

BMVert *BM_vert_at_index(BMesh *bm, const int index)
{
	BLI_assert((index >= 0) && (index < bm->totvert));
	BLI_assert((bm->elem_table_dirty & BM_VERT) == 0);
	return bm->vtable[index];
}

BMEdge *BM_edge_at_index(BMesh *bm, const int index)
{
	BLI_assert((index >= 0) && (index < bm->totedge));
	BLI_assert((bm->elem_table_dirty & BM_EDGE) == 0);
	return bm->etable[index];
}

BMFace *BM_face_at_index(BMesh *bm, const int index)
{
	BLI_assert((index >= 0) && (index < bm->totface));
	BLI_assert((bm->elem_table_dirty & BM_FACE) == 0);
	return bm->ftable[index];
}


BMVert *BM_vert_at_index_find(BMesh *bm, const int index)
{
	return BLI_mempool_findelem(bm->vpool, index);
}

BMEdge *BM_edge_at_index_find(BMesh *bm, const int index)
{
	return BLI_mempool_findelem(bm->epool, index);
}

BMFace *BM_face_at_index_find(BMesh *bm, const int index)
{
	return BLI_mempool_findelem(bm->fpool, index);
}


/**
 * Return the amount of element of type 'type' in a given bmesh.
 */
int BM_mesh_elem_count(BMesh *bm, const char htype)
{
	BLI_assert((htype & ~BM_ALL_NOLOOP) == 0);

	switch (htype) {
		case BM_VERT: return bm->totvert;
		case BM_EDGE: return bm->totedge;
		case BM_FACE: return bm->totface;
		default:
		{
			BLI_assert(0);
			return 0;
		}
	}
}

/**
 * Remaps the vertices, edges and/or faces of the bmesh as indicated by vert/edge/face_idx arrays
 * (xxx_idx[org_index] = new_index).
 *
 * A NULL array means no changes.
 *
 * Note: - Does not mess with indices, just sets elem_index_dirty flag.
 *       - For verts/edges/faces only (as loops must remain "ordered" and "aligned"
 *         on a per-face basis...).
 *
 * WARNING: Be careful if you keep pointers to affected BM elements, or arrays, when using this func!
 */
void BM_mesh_remap(BMesh *bm, unsigned int *vert_idx, unsigned int *edge_idx, unsigned int *face_idx)
{
	/* Mapping old to new pointers. */
	GHash *vptr_map = NULL, *eptr_map = NULL, *fptr_map = NULL;
	BMIter iter, iterl;
	BMVert *ve;
	BMEdge *ed;
	BMFace *fa;
	BMLoop *lo;

	if (!(vert_idx || edge_idx || face_idx))
		return;

	/* Remap Verts */
	if (vert_idx) {
		BMVert **verts_pool, *verts_copy, **vep;
		int i, totvert = bm->totvert;
		unsigned int *new_idx = NULL;

		/* Init the old-to-new vert pointers mapping */
		vptr_map = BLI_ghash_ptr_new_ex("BM_mesh_remap vert pointers mapping", bm->totvert);

		/* Make a copy of all vertices. */
		verts_pool = MEM_callocN(sizeof(BMVert *) * totvert, "BM_mesh_remap verts pool");
		BM_iter_as_array(bm, BM_VERTS_OF_MESH, NULL, (void **)verts_pool, totvert);
		verts_copy = MEM_mallocN(sizeof(BMVert) * totvert, "BM_mesh_remap verts copy");
		for (i = totvert, ve = verts_copy + totvert - 1, vep = verts_pool + totvert - 1; i--; ve--, vep--) {
			*ve = **vep;
/*			printf("*vep: %p, verts_pool[%d]: %p\n", *vep, i, verts_pool[i]);*/
		}

		/* Copy back verts to their new place, and update old2new pointers mapping. */
		new_idx = vert_idx + totvert - 1;
		ve = verts_copy + totvert - 1;
		vep = verts_pool + totvert - 1; /* old, org pointer */
		for (i = totvert; i--; new_idx--, ve--, vep--) {
			BMVert *new_vep = verts_pool[*new_idx];
			*new_vep = *ve;
/*			printf("mapping vert from %d to %d (%p/%p to %p)\n", i, *new_idx, *vep, verts_pool[i], new_vep);*/
			BLI_ghash_insert(vptr_map, (void *)*vep, (void *)new_vep);
		}
		bm->elem_index_dirty |= BM_VERT;

		MEM_freeN(verts_pool);
		MEM_freeN(verts_copy);
	}

	/* Remap Edges */
	if (edge_idx) {
		BMEdge **edges_pool, *edges_copy, **edp;
		int i, totedge = bm->totedge;
		unsigned int *new_idx = NULL;

		/* Init the old-to-new vert pointers mapping */
		eptr_map = BLI_ghash_ptr_new_ex("BM_mesh_remap edge pointers mapping", bm->totedge);

		/* Make a copy of all vertices. */
		edges_pool = MEM_callocN(sizeof(BMEdge *) * totedge, "BM_mesh_remap edges pool");
		BM_iter_as_array(bm, BM_EDGES_OF_MESH, NULL, (void **)edges_pool, totedge);
		edges_copy = MEM_mallocN(sizeof(BMEdge) * totedge, "BM_mesh_remap edges copy");
		for (i = totedge, ed = edges_copy + totedge - 1, edp = edges_pool + totedge - 1; i--; ed--, edp--) {
			*ed = **edp;
		}

		/* Copy back verts to their new place, and update old2new pointers mapping. */
		new_idx = edge_idx + totedge - 1;
		ed = edges_copy + totedge - 1;
		edp = edges_pool + totedge - 1; /* old, org pointer */
		for (i = totedge; i--; new_idx--, ed--, edp--) {
			BMEdge *new_edp = edges_pool[*new_idx];
			*new_edp = *ed;
			BLI_ghash_insert(eptr_map, (void *)*edp, (void *)new_edp);
/*			printf("mapping edge from %d to %d (%p/%p to %p)\n", i, *new_idx, *edp, edges_pool[i], new_edp);*/
		}
		bm->elem_index_dirty |= BM_EDGE;

		MEM_freeN(edges_pool);
		MEM_freeN(edges_copy);
	}

	/* Remap Faces */
	if (face_idx) {
		BMFace **faces_pool, *faces_copy, **fap;
		int i, totface = bm->totface;
		unsigned int *new_idx = NULL;

		/* Init the old-to-new vert pointers mapping */
		fptr_map = BLI_ghash_ptr_new_ex("BM_mesh_remap face pointers mapping", bm->totface);

		/* Make a copy of all vertices. */
		faces_pool = MEM_callocN(sizeof(BMFace *) * totface, "BM_mesh_remap faces pool");
		BM_iter_as_array(bm, BM_FACES_OF_MESH, NULL, (void **)faces_pool, totface);
		faces_copy = MEM_mallocN(sizeof(BMFace) * totface, "BM_mesh_remap faces copy");
		for (i = totface, fa = faces_copy + totface - 1, fap = faces_pool + totface - 1; i--; fa--, fap--) {
			*fa = **fap;
		}

		/* Copy back verts to their new place, and update old2new pointers mapping. */
		new_idx = face_idx + totface - 1;
		fa = faces_copy + totface - 1;
		fap = faces_pool + totface - 1; /* old, org pointer */
		for (i = totface; i--; new_idx--, fa--, fap--) {
			BMFace *new_fap = faces_pool[*new_idx];
			*new_fap = *fa;
			BLI_ghash_insert(fptr_map, (void *)*fap, (void *)new_fap);
		}

		bm->elem_index_dirty |= BM_FACE;

		MEM_freeN(faces_pool);
		MEM_freeN(faces_copy);
	}

	/* And now, fix all vertices/edges/faces/loops pointers! */
	/* Verts' pointers, only edge pointers... */
	if (eptr_map) {
		BM_ITER_MESH (ve, &iter, bm, BM_VERTS_OF_MESH) {
/*			printf("Vert e: %p -> %p\n", ve->e, BLI_ghash_lookup(eptr_map, (const void *)ve->e));*/
			ve->e = BLI_ghash_lookup(eptr_map, (const void *)ve->e);
		}
	}

	/* Edges' pointers, only vert pointers (as we don't mess with loops!), and - ack! - edge pointers,
	 * as we have to handle disklinks... */
	if (vptr_map || eptr_map) {
		BM_ITER_MESH (ed, &iter, bm, BM_EDGES_OF_MESH) {
			if (vptr_map) {
/*				printf("Edge v1: %p -> %p\n", ed->v1, BLI_ghash_lookup(vptr_map, (const void *)ed->v1));*/
/*				printf("Edge v2: %p -> %p\n", ed->v2, BLI_ghash_lookup(vptr_map, (const void *)ed->v2));*/
				ed->v1 = BLI_ghash_lookup(vptr_map, (const void *)ed->v1);
				ed->v2 = BLI_ghash_lookup(vptr_map, (const void *)ed->v2);
			}
			if (eptr_map) {
/*				printf("Edge v1_disk_link prev: %p -> %p\n", ed->v1_disk_link.prev,*/
/*				       BLI_ghash_lookup(eptr_map, (const void *)ed->v1_disk_link.prev));*/
/*				printf("Edge v1_disk_link next: %p -> %p\n", ed->v1_disk_link.next,*/
/*				       BLI_ghash_lookup(eptr_map, (const void *)ed->v1_disk_link.next));*/
/*				printf("Edge v2_disk_link prev: %p -> %p\n", ed->v2_disk_link.prev,*/
/*				       BLI_ghash_lookup(eptr_map, (const void *)ed->v2_disk_link.prev));*/
/*				printf("Edge v2_disk_link next: %p -> %p\n", ed->v2_disk_link.next,*/
/*				       BLI_ghash_lookup(eptr_map, (const void *)ed->v2_disk_link.next));*/
				ed->v1_disk_link.prev = BLI_ghash_lookup(eptr_map, (const void *)ed->v1_disk_link.prev);
				ed->v1_disk_link.next = BLI_ghash_lookup(eptr_map, (const void *)ed->v1_disk_link.next);
				ed->v2_disk_link.prev = BLI_ghash_lookup(eptr_map, (const void *)ed->v2_disk_link.prev);
				ed->v2_disk_link.next = BLI_ghash_lookup(eptr_map, (const void *)ed->v2_disk_link.next);
			}
		}
	}

	/* Faces' pointers (loops, in fact), always needed... */
	BM_ITER_MESH (fa, &iter, bm, BM_FACES_OF_MESH) {
		BM_ITER_ELEM (lo, &iterl, fa, BM_LOOPS_OF_FACE) {
			if (vptr_map) {
/*				printf("Loop v: %p -> %p\n", lo->v, BLI_ghash_lookup(vptr_map, (const void *)lo->v));*/
				lo->v = BLI_ghash_lookup(vptr_map, (const void *)lo->v);
			}
			if (eptr_map) {
/*				printf("Loop e: %p -> %p\n", lo->e, BLI_ghash_lookup(eptr_map, (const void *)lo->e));*/
				lo->e = BLI_ghash_lookup(eptr_map, (const void *)lo->e);
			}
			if (fptr_map) {
/*				printf("Loop f: %p -> %p\n", lo->f, BLI_ghash_lookup(fptr_map, (const void *)lo->f));*/
				lo->f = BLI_ghash_lookup(fptr_map, (const void *)lo->f);
			}
		}
	}

	if (vptr_map)
		BLI_ghash_free(vptr_map, NULL, NULL);
	if (eptr_map)
		BLI_ghash_free(eptr_map, NULL, NULL);
	if (fptr_map)
		BLI_ghash_free(fptr_map, NULL, NULL);
}

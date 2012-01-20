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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) Blender Foundation, 2002-2009
 * All rights reserved.
 *
 * Contributor(s): Antony Riakiotakis
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 * UV Sculpt tools
 *
 */

/** \file blender/editors/sculpt_paint/sculpt_uv.c
 *  \ingroup edsculpt
 */


#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"
#include "BLI_editVert.h"
#include "BLI_math.h"
#include "BLI_ghash.h"

#include "DNA_object_types.h"
#include "DNA_scene_types.h"
#include "DNA_brush_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_brush.h"
#include "BKE_paint.h"
#include "BKE_context.h"
#include "BKE_main.h"
#include "BKE_depsgraph.h"
#include "BKE_mesh.h"
#include "BKE_customdata.h"

#include "ED_screen.h"
#include "ED_image.h"
#include "ED_mesh.h"

#include "WM_api.h"
#include "WM_types.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "paint_intern.h"
#include "uvedit_intern.h"

#include "UI_view2d.h"

#define MARK_BOUNDARY	1

typedef struct UvAdjacencyElement {
	/* pointer to original uvelement */
	UvElement *element;
	/* uv pointer for convenience. Caution, this points to the original UVs! */
	float *uv;
	/* general use flag (Used to check if Element is boundary here) */
	char flag;
} UvAdjacencyElement;

typedef struct UvEdge {
	unsigned int uv1;
	unsigned int uv2;
	/* general use flag (Used to check if edge is boundary here, and propagates to adjacency elements) */
	char flag;
}UvEdge;

typedef struct UVInitialStrokeElement{
	/* index to unique uv */
	int uv;

	/* strength of brush on initial position */
	float strength;

	/* initial uv position */
	float initial_uv[2];
}UVInitialStrokeElement;

typedef struct UVInitialStroke{
	/* Initial Selection,for grab brushes for instance */
	UVInitialStrokeElement *initialSelection;

	/* total initially selected UVs*/
	int totalInitialSelected;

	/* initial mouse coordinates */
	float init_coord[2];
}UVInitialStroke;


/* custom data for uv smoothing brush */
typedef struct UvSculptData{
	/* Contains the first of each set of coincident uvs.
	 * These will be used to perform smoothing on and propagate the changes
	 * to their coincident uvs */
	UvAdjacencyElement *uv;

	/* ...Is what it says */
	int totalUniqueUvs;

	/* Edges used for adjacency info, used with laplacian smoothing */
	UvEdge *uvedges;

	/* Need I say more? */
	int totalUvEdges;

	/* data for initial stroke, used by tools like grab */
	UVInitialStroke *initial_stroke;

	/* Timer to be used for airbrush-type brush */
	wmTimer *timer;

	/* To determine quickly adjacent uvs */
	UvElementMap *elementMap;

	/* uvsmooth Paint for fast reference */
	Paint *uvsculpt;
}UvSculptData;

/*********** Improved Laplacian Relaxation Operator ************************/
/* Original code by Raul Fernandez Hernandez "farsthary"                   *
 * adapted to uv smoothing by Antony Riakiatakis                           *
 ***************************************************************************/

typedef struct Temp_UvData{
	float sum_co[2], p[2], b[2], sum_b[2];
	int ncounter;
}Temp_UVData;



void HC_relaxation_iteration_uv(EditMesh *em, UvSculptData *sculptdata, float mouse_coord[2], float alpha, float radius, float aspectRatio){
	Temp_UVData *tmp_uvdata;
	float diff[2];
	int i;
	float radius_root = sqrt(radius);
	Brush *brush = paint_brush(sculptdata->uvsculpt);

	tmp_uvdata = (Temp_UVData *)MEM_callocN(sculptdata->totalUniqueUvs * sizeof(Temp_UVData), "Temporal data");

	/* counting neighbors */
	for (i = 0; i < sculptdata->totalUvEdges; i++){
		UvEdge *tmpedge = sculptdata->uvedges+i;
		tmp_uvdata[tmpedge->uv1].ncounter++;
		tmp_uvdata[tmpedge->uv2].ncounter++;

		add_v2_v2(tmp_uvdata[tmpedge->uv2].sum_co, sculptdata->uv[tmpedge->uv1].uv);
		add_v2_v2(tmp_uvdata[tmpedge->uv1].sum_co, sculptdata->uv[tmpedge->uv2].uv);
	}

	for (i = 0; i < sculptdata->totalUniqueUvs; i++){
		copy_v2_v2(diff,tmp_uvdata[i].sum_co);
		mul_v2_fl(diff,1.f/tmp_uvdata[i].ncounter);
		copy_v2_v2(tmp_uvdata[i].p,diff);

		tmp_uvdata[i].b[0] = diff[0] - sculptdata->uv[i].uv[0];
		tmp_uvdata[i].b[1] = diff[1] - sculptdata->uv[i].uv[1];
	}

	for (i = 0; i < sculptdata->totalUvEdges; i++){
		UvEdge *tmpedge = sculptdata->uvedges+i;
		add_v2_v2(tmp_uvdata[tmpedge->uv1].sum_b, tmp_uvdata[tmpedge->uv2].b);
		add_v2_v2(tmp_uvdata[tmpedge->uv2].sum_b, tmp_uvdata[tmpedge->uv1].b);
	}

	for (i = 0; i < sculptdata->totalUniqueUvs; i++){
		float dist;
		/* This is supposed to happen only if "Pin Edges" is on, since we have initialization on stroke start
		 * If ever uv brushes get their own mode we should check for toolsettings option too */
		if((sculptdata->uv[i].flag & MARK_BOUNDARY)){
			continue;
		}

		sub_v2_v2v2(diff, sculptdata->uv[i].uv, mouse_coord);
		diff[1] /= aspectRatio;
		if((dist = dot_v2v2(diff, diff)) <= radius){
			UvElement *element;
			float strength;
			strength = alpha*brush_curve_strength(brush, sqrt(dist), radius_root);

			sculptdata->uv[i].uv[0] = (1.0-strength)*sculptdata->uv[i].uv[0] + strength*(tmp_uvdata[i].p[0] - 0.5f*(tmp_uvdata[i].b[0] + tmp_uvdata[i].sum_b[0]/tmp_uvdata[i].ncounter));
			sculptdata->uv[i].uv[1] = (1.0-strength)*sculptdata->uv[i].uv[1] + strength*(tmp_uvdata[i].p[1] - 0.5f*(tmp_uvdata[i].b[1] + tmp_uvdata[i].sum_b[1]/tmp_uvdata[i].ncounter));

			for(element = sculptdata->uv[i].element; element; element = element->next){
				MTFace *mt;
				if(element->separate && element != sculptdata->uv[i].element)
					break;
				mt = CustomData_em_get(&em->fdata, element->face->data, CD_MTFACE);
				copy_v2_v2(mt->uv[element->tfindex], sculptdata->uv[i].uv);
			}
		}
	}

	MEM_freeN(tmp_uvdata);

	return;
}

static void laplacian_relaxation_iteration_uv(EditMesh *em, UvSculptData *sculptdata, float mouse_coord[2], float alpha, float radius, float aspectRatio)
{
	Temp_UVData *tmp_uvdata;
	float diff[2];
	int i;
	float radius_root = sqrt(radius);
	Brush *brush = paint_brush(sculptdata->uvsculpt);

	tmp_uvdata = (Temp_UVData *)MEM_callocN(sculptdata->totalUniqueUvs * sizeof(Temp_UVData), "Temporal data");

	/* counting neighbors */
	for (i = 0; i < sculptdata->totalUvEdges; i++){
		UvEdge *tmpedge = sculptdata->uvedges+i;
		tmp_uvdata[tmpedge->uv1].ncounter++;
		tmp_uvdata[tmpedge->uv2].ncounter++;

		add_v2_v2(tmp_uvdata[tmpedge->uv2].sum_co, sculptdata->uv[tmpedge->uv1].uv);
		add_v2_v2(tmp_uvdata[tmpedge->uv1].sum_co, sculptdata->uv[tmpedge->uv2].uv);
	}

	/* Original Lacplacian algorithm included removal of normal component of translation. here it is not
	 * needed since we translate along the UV plane always.*/
	for (i = 0; i < sculptdata->totalUniqueUvs; i++){
		copy_v2_v2(tmp_uvdata[i].p, tmp_uvdata[i].sum_co);
		mul_v2_fl(tmp_uvdata[i].p, 1.f/tmp_uvdata[i].ncounter);
	}

	for (i = 0; i < sculptdata->totalUniqueUvs; i++){
		float dist;
		/* This is supposed to happen only if "Pin Edges" is on, since we have initialization on stroke start
		 * If ever uv brushes get their own mode we should check for toolsettings option too */
		if((sculptdata->uv[i].flag & MARK_BOUNDARY)){
			continue;
		}

		sub_v2_v2v2(diff, sculptdata->uv[i].uv, mouse_coord);
		diff[1] /= aspectRatio;
		if((dist = dot_v2v2(diff, diff)) <= radius){
			UvElement *element;
			float strength;
			strength = alpha*brush_curve_strength(brush, sqrt(dist), radius_root);

			sculptdata->uv[i].uv[0] = (1.0-strength)*sculptdata->uv[i].uv[0] + strength*tmp_uvdata[i].p[0];
			sculptdata->uv[i].uv[1] = (1.0-strength)*sculptdata->uv[i].uv[1] + strength*tmp_uvdata[i].p[1];

			for(element = sculptdata->uv[i].element; element; element = element->next){
				MTFace *mt;
				if(element->separate && element != sculptdata->uv[i].element)
					break;
				mt = CustomData_em_get(&em->fdata, element->face->data, CD_MTFACE);
				copy_v2_v2(mt->uv[element->tfindex], sculptdata->uv[i].uv);
			}
		}
	}

	MEM_freeN(tmp_uvdata);

	return;
}


static void uv_sculpt_stroke_apply(bContext *C, wmOperator *op, wmEvent *event, Object *obedit)
{
	float co[2], radius, radius_root;
	Scene *scene = CTX_data_scene(C);
	ARegion *ar = CTX_wm_region(C);
	EditMesh *em = BKE_mesh_get_editmesh(obedit->data);
	unsigned int tool;
	UvSculptData *sculptdata = (UvSculptData *)op->customdata;
	SpaceImage *sima;
	int invert;
	int width, height;
	float aspectRatio;
	float alpha, zoomx, zoomy;
	Brush *brush = paint_brush(sculptdata->uvsculpt);
	ToolSettings *toolsettings = CTX_data_tool_settings(C);
	tool = RNA_boolean_get(op->ptr, "temp_relax")? UV_SCULPT_TOOL_RELAX : toolsettings->uv_sculpt_tool;

	invert = RNA_boolean_get(op->ptr, "invert")? -1 : 1;
	alpha = brush_alpha(scene, brush);
	UI_view2d_region_to_view(&ar->v2d, event->mval[0], event->mval[1], &co[0], &co[1]);

	sima = CTX_wm_space_image(C);
	ED_space_image_size(sima, &width, &height);
	ED_space_image_zoom(sima, ar, &zoomx, &zoomy);

	radius = brush_size(scene, brush)/(width*zoomx);
	aspectRatio = width/(float)height;

	/* We will compare squares to save some computation */
	radius = radius*radius;
	radius_root = sqrt(radius);

	/*
	 * Pinch Tool
	 */
	if(tool == UV_SCULPT_TOOL_PINCH){
		int i;
		alpha *= invert;
		for (i = 0; i < sculptdata->totalUniqueUvs; i++){
			float dist, diff[2];
			/* This is supposed to happen only if "Lock Borders" is on, since we have initialization on stroke start
			 * If ever uv brushes get their own mode we should check for toolsettings option too */
			if(sculptdata->uv[i].flag & MARK_BOUNDARY){
				continue;
			}

			sub_v2_v2v2(diff, sculptdata->uv[i].uv, co);
			diff[1] /= aspectRatio;
			if((dist = dot_v2v2(diff, diff)) <= radius){
				UvElement *element;
				float strength;
				strength = alpha*brush_curve_strength(brush, sqrt(dist), radius_root);
				normalize_v2(diff);

				sculptdata->uv[i].uv[0] -= strength*diff[0]*0.001;
				sculptdata->uv[i].uv[1] -= strength*diff[1]*0.001;

				for(element = sculptdata->uv[i].element; element; element = element->next){
					MTFace *mt;
					if(element->separate && element != sculptdata->uv[i].element)
						break;
					mt = CustomData_em_get(&em->fdata, element->face->data, CD_MTFACE);
					copy_v2_v2(mt->uv[element->tfindex], sculptdata->uv[i].uv);
				}
			}
		}
	}

	/*
	 * Smooth Tool
	 */
	else if(tool == UV_SCULPT_TOOL_RELAX){
		unsigned int method = toolsettings->uv_relax_method;
		if(method == UV_SCULPT_TOOL_RELAX_HC){
			HC_relaxation_iteration_uv(em, sculptdata, co, alpha, radius, aspectRatio);
		}else{
			laplacian_relaxation_iteration_uv(em, sculptdata, co, alpha, radius, aspectRatio);
		}
	}

	/*
	 * Grab Tool
	 */
	else if(tool == UV_SCULPT_TOOL_GRAB){
		int i;
		float diff[2];
		sub_v2_v2v2(diff, co, sculptdata->initial_stroke->init_coord);

		for(i = 0; i < sculptdata->initial_stroke->totalInitialSelected; i++ ){
			UvElement *element;
			int uvindex = sculptdata->initial_stroke->initialSelection[i].uv;
			float strength = sculptdata->initial_stroke->initialSelection[i].strength;
			sculptdata->uv[uvindex].uv[0] = sculptdata->initial_stroke->initialSelection[i].initial_uv[0] + strength*diff[0];
			sculptdata->uv[uvindex].uv[1] = sculptdata->initial_stroke->initialSelection[i].initial_uv[1] + strength*diff[1];

			for(element = sculptdata->uv[uvindex].element; element; element = element->next){
				MTFace *mt;
				if(element->separate && element != sculptdata->uv[uvindex].element)
					break;
				mt = CustomData_em_get(&em->fdata, element->face->data, CD_MTFACE);
				copy_v2_v2(mt->uv[element->tfindex], sculptdata->uv[uvindex].uv);
			}
		}
	}

	BKE_mesh_end_editmesh(obedit->data, em);
}


static void uv_sculpt_stroke_exit(bContext *C, wmOperator *op)
{
	UvSculptData *data = op->customdata;
	if(data->timer){
		WM_event_remove_timer(CTX_wm_manager(C), CTX_wm_window(C), data->timer);
	}
	if(data->elementMap)
	{
		EM_free_uv_element_map(data->elementMap);
	}
	if(data->uv){
		MEM_freeN(data->uv);
	}
	if(data->uvedges){
		MEM_freeN(data->uvedges);
	}
	if(data->initial_stroke){
		if(data->initial_stroke->initialSelection){
			MEM_freeN(data->initial_stroke->initialSelection);
		}
		MEM_freeN(data->initial_stroke);
	}

	MEM_freeN(data);
	op->customdata = NULL;
}

static int get_uv_element_offset_from_face(UvElementMap *map, EditFace *efa, int index, int island_index, int doIslands){
	UvElement *element = ED_get_uv_element(map, efa, index);
	if(!element || (doIslands && element->island != island_index)){
		return -1;
	}
	return element - map->buf;
}


static unsigned int	uv_edge_hash(const void *key){
	UvEdge *edge = (UvEdge *)key;
	return 
		BLI_ghashutil_inthash(SET_INT_IN_POINTER(edge->uv2)) +
		BLI_ghashutil_inthash(SET_INT_IN_POINTER(edge->uv1));
}

static int uv_edge_compare(const void *a, const void *b){
	UvEdge *edge1 = (UvEdge *)a;
	UvEdge *edge2 = (UvEdge *)b;

	if((edge1->uv1 == edge2->uv1) && (edge1->uv2 == edge2->uv2)){
		return 0;
	}
	return 1;
}


static UvSculptData *uv_sculpt_stroke_init(bContext *C, wmOperator *op, wmEvent *event)
{
	Scene *scene = CTX_data_scene(C);
	Object *obedit = CTX_data_edit_object(C);
	ToolSettings *ts = scene->toolsettings;
	UvSculptData *data = MEM_callocN(sizeof(*data), "UV Smooth Brush Data");
	EditMesh *em = BKE_mesh_get_editmesh(obedit->data);

	op->customdata = data;

	if(data){
		int counter = 0, i;
		ARegion *ar= CTX_wm_region(C);
		float co[2];
		EditFace *efa;
		UvEdge *edges;
		GHash *edgeHash;
		GHashIterator* ghi;
		MTFace *mt;
		int do_island_optimization = !(ts->uv_sculpt_settings & UV_SCULPT_ALL_ISLANDS);
		int island_index = 0;
		/* Holds, for each UvElement in elementMap, a pointer to its unique uv.*/
		int *uniqueUv;

		data->uvsculpt = &ts->uvsculpt->paint;

		if(do_island_optimization){
			/* We will need island information */
			if(ts->uv_flag & UV_SYNC_SELECTION){
				data->elementMap = EM_make_uv_element_map(em, 0, 1);
			}else{
				data->elementMap = EM_make_uv_element_map(em, 1, 1);
			}
		}else {
			if(ts->uv_flag & UV_SYNC_SELECTION){
				data->elementMap = EM_make_uv_element_map(em, 0, 0);
			}else{
				data->elementMap = EM_make_uv_element_map(em, 1, 0);
			}
		}

		if(!data->elementMap){
			uv_sculpt_stroke_exit(C, op);
			return NULL;
		}

		/* Mouse coordinates, useful for some functions like grab and sculpt all islands */
		UI_view2d_region_to_view(&ar->v2d, event->mval[0], event->mval[1], &co[0], &co[1]);

		/* we need to find the active island here */
		if(do_island_optimization){
			UvElement *element;
			NearestHit hit;
			Image *ima= CTX_data_edit_image(C);
			uv_find_nearest_vert(scene, ima, em, co, NULL, &hit);

			element = ED_get_uv_element(data->elementMap, hit.efa, hit.uv);
			island_index = element->island;
		}


		/* Count 'unique' uvs */
		for(i = 0; i < data->elementMap->totalUVs; i++){
			if(data->elementMap->buf[i].separate
			&& (!do_island_optimization || data->elementMap->buf[i].island == island_index)){
				counter++;
			}
		}

		/* Allocate the unique uv buffers */
		data->uv = MEM_mallocN(sizeof(*data->uv)*counter, "uv_brush_unique_uvs");
		uniqueUv = MEM_mallocN(sizeof(*uniqueUv)*data->elementMap->totalUVs, "uv_brush_unique_uv_map");
		edgeHash = BLI_ghash_new(uv_edge_hash, uv_edge_compare, "uv_brush_edge_hash");
		/* we have at most totalUVs edges */
		edges = MEM_mallocN(sizeof(*edges)*data->elementMap->totalUVs, "uv_brush_all_edges");
		if(!data->uv || !uniqueUv || !edgeHash || !edges){
			if(edges){
				MEM_freeN(edges);
			}
			if(uniqueUv){
				MEM_freeN(uniqueUv);
			}
			if(edgeHash){
				MEM_freeN(edgeHash);
			}
			uv_sculpt_stroke_exit(C, op);
			return NULL;
		}

		data->totalUniqueUvs = counter;
		/* So that we can use this as index for the UvElements */
		counter = -1;
		/* initialize the unique UVs */
		for(i = 0; i < em->totvert; i++){
			UvElement *element = data->elementMap->vert[i];
			for(; element; element = element->next){
				if(element->separate){
					if(do_island_optimization && (element->island != island_index)){
						/* skip this uv if not on the active island */
						for(; element->next && !(element->next->separate); element = element->next)
							;
						continue;
					}
					efa = element->face;
					mt = CustomData_em_get(&em->fdata, efa->data, CD_MTFACE);

					counter++;
					data->uv[counter].element = element;
					data->uv[counter].flag = 0;
					data->uv[counter].uv = mt->uv[element->tfindex];
				}
				/* pointer arithmetic to the rescue, as always :)*/
				uniqueUv[element - data->elementMap->buf] = counter;
			}
		}

		/* Now, on to generate our uv connectivity data */
		for(efa = em->faces.first, counter = 0; efa; efa = efa->next){
			int nverts = efa->v4 ? 4 : 3;
			for(i = 0; i < nverts; i++){
				int offset1, itmp1 = get_uv_element_offset_from_face(data->elementMap, efa, i, island_index, do_island_optimization);
				int offset2, itmp2 = get_uv_element_offset_from_face(data->elementMap, efa, (i+1)%nverts, island_index, do_island_optimization);

				/* Skip edge if not found(unlikely) or not on valid island */
				if(itmp1 == -1 || itmp2 == -1)
					continue;

				offset1 = uniqueUv[itmp1];
				offset2 = uniqueUv[itmp2];

				edges[counter].flag = 0;
				/* using an order policy, sort uvs according to address space. This avoids
				 * Having two different UvEdges with the same uvs on different positions  */
				if(offset1 < offset2){
					edges[counter].uv1 = offset1;
					edges[counter].uv2 = offset2;
				}
				else{
					edges[counter].uv1 = offset2;
					edges[counter].uv2 = offset1;
				}
				/* Hack! Set the value of the key to its flag. Now we can set the flag when an edge exists twice :) */
				if(BLI_ghash_haskey(edgeHash, &edges[counter])){
					char *flag = BLI_ghash_lookup(edgeHash, &edges[counter]);
					*flag = 1;
				}
				else{
					/* Hack mentioned */
					BLI_ghash_insert(edgeHash, &edges[counter], &edges[counter].flag);
				}
				counter++;
			}
		}
		MEM_freeN(uniqueUv);

		/* Allocate connectivity data, we allocate edges once */
		data->uvedges = MEM_mallocN(sizeof(*data->uvedges)*BLI_ghash_size(edgeHash), "uv_brush_edge_connectivity_data");
		if(!data->uvedges){
			BLI_ghash_free(edgeHash, NULL, NULL);
			MEM_freeN(edges);
			uv_sculpt_stroke_exit(C, op);
			return NULL;
		}
		ghi = BLI_ghashIterator_new(edgeHash);
		if(!ghi){
			BLI_ghash_free(edgeHash, NULL, NULL);
			MEM_freeN(edges);
			uv_sculpt_stroke_exit(C, op);
			return NULL;
		}
		/* fill the edges with data */
		for(i = 0; !BLI_ghashIterator_isDone(ghi); BLI_ghashIterator_step(ghi)){
			data->uvedges[i++] = *((UvEdge *)BLI_ghashIterator_getKey(ghi));
		}
		data->totalUvEdges = BLI_ghash_size(edgeHash);

		/* cleanup temporary stuff */
		BLI_ghashIterator_free(ghi);
		BLI_ghash_free(edgeHash, NULL, NULL);
		MEM_freeN(edges);

		/* transfer boundary edge property to uvs */
		if(ts->uv_sculpt_settings & UV_SCULPT_LOCK_BORDERS){
			for(i = 0; i < data->totalUvEdges; i++){
				if(!data->uvedges[i].flag){
					data->uv[data->uvedges[i].uv1].flag |= MARK_BOUNDARY;
					data->uv[data->uvedges[i].uv2].flag |= MARK_BOUNDARY;
				}
			}
		}

		/* Allocate initial selection for grab tool */
		if(ts->uv_sculpt_tool == UV_SCULPT_TOOL_GRAB){
			float radius, radius_root;
			UvSculptData *sculptdata = (UvSculptData *)op->customdata;
			SpaceImage *sima;
			int width, height;
			float aspectRatio;
			float alpha, zoomx, zoomy;
			Brush *brush = paint_brush(sculptdata->uvsculpt);

			alpha = brush_alpha(scene, brush);

			radius = brush_size(scene, brush);
			sima = CTX_wm_space_image(C);
			ED_space_image_size(sima, &width, &height);
			ED_space_image_zoom(sima, ar, &zoomx, &zoomy);

			aspectRatio = width/(float)height;
			radius /= (width*zoomx);
			radius = radius*radius;
			radius_root = sqrt(radius);

			/* Allocate selection stack */
			data->initial_stroke = MEM_mallocN(sizeof(*data->initial_stroke), "uv_sculpt_initial_stroke");
			if(!data->initial_stroke){
				uv_sculpt_stroke_exit(C, op);
			}
			data->initial_stroke->initialSelection = MEM_mallocN(sizeof(*data->initial_stroke->initialSelection)*data->totalUniqueUvs, "uv_sculpt_initial_selection");
			if(!data->initial_stroke->initialSelection){
				uv_sculpt_stroke_exit(C, op);
			}

			copy_v2_v2(data->initial_stroke->init_coord, co);

			counter = 0;

			for(i = 0; i < data->totalUniqueUvs; i++){
				float dist, diff[2];
				if(data->uv[i].flag & MARK_BOUNDARY){
					continue;
				}

				sub_v2_v2v2(diff, data->uv[i].uv, co);
				diff[1] /= aspectRatio;
				if((dist = dot_v2v2(diff, diff)) <= radius){
					float strength;
					strength = alpha*brush_curve_strength(brush, sqrt(dist), radius_root);

					data->initial_stroke->initialSelection[counter].uv = i;
					data->initial_stroke->initialSelection[counter].strength = strength;
					copy_v2_v2(data->initial_stroke->initialSelection[counter].initial_uv, data->uv[i].uv);
					counter++;
				}
			}

			data->initial_stroke->totalInitialSelected = counter;
		}
	}

	return op->customdata;
}

static int uv_sculpt_stroke_invoke(bContext *C, wmOperator *op, wmEvent *event)
{
	UvSculptData *data;
	Object *obedit = CTX_data_edit_object(C);

	if(!(data = uv_sculpt_stroke_init(C, op, event))) {
		return OPERATOR_CANCELLED;
	}

	uv_sculpt_stroke_apply(C, op, event, obedit);

	data->timer= WM_event_add_timer(CTX_wm_manager(C), CTX_wm_window(C), TIMER, 0.001f);

	if(!data->timer){
		uv_sculpt_stroke_exit(C, op);
		return OPERATOR_CANCELLED;
	}
	WM_event_add_modal_handler(C, op);

	return OPERATOR_RUNNING_MODAL;
}


static int uv_sculpt_stroke_modal(bContext *C, wmOperator *op, wmEvent *event)
{
	UvSculptData *data = (UvSculptData *)op->customdata;
	Object *obedit = CTX_data_edit_object(C);

	switch(event->type) {
		case LEFTMOUSE:
		case MIDDLEMOUSE:
		case RIGHTMOUSE:
			uv_sculpt_stroke_exit(C, op);
			return OPERATOR_FINISHED;

		case MOUSEMOVE:
		case INBETWEEN_MOUSEMOVE:
			uv_sculpt_stroke_apply(C, op, event, obedit);
			break;
		case TIMER:
			if(event->customdata == data->timer)
				uv_sculpt_stroke_apply(C, op, event, obedit);
			break;
		default:
			return OPERATOR_RUNNING_MODAL;
	}

	ED_region_tag_redraw(CTX_wm_region(C));
	WM_event_add_notifier(C, NC_GEOM|ND_DATA, obedit->data);
	DAG_id_tag_update(obedit->data, 0);
	return OPERATOR_RUNNING_MODAL;
}

void SCULPT_OT_uv_sculpt_stroke(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Sculpt UVs";
	ot->description = "Sculpt UVs using a brush";
	ot->idname = "SCULPT_OT_uv_sculpt_stroke";

	/* api callbacks */
	ot->invoke = uv_sculpt_stroke_invoke;
	ot->modal = uv_sculpt_stroke_modal;
	ot->poll = uv_sculpt_poll;

	/* flags */
	ot->flag = OPTYPE_REGISTER|OPTYPE_UNDO;

	/* props */
	RNA_def_boolean(ot->srna, "invert", 0, "Invert", "Inverts the operator");
	RNA_def_boolean(ot->srna, "temp_relax", 0, "Relax", "Relax Tool");
}

/*
 * $Id: ama.c 36773 2011-08-13 13:46:00Z ruesp83 $
 *
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
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Fabio Russo
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/ama.c
 *  \ingroup bke
 */

#include "MEM_guardedalloc.h"
#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_rand.h"
#include "BLI_listbase.h"
#include "BLI_edgehash.h"
#include "DNA_group_types.h"
#include "DNA_object_types.h"
#include "DNA_curve_types.h"
#include "DNA_modifier_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_scene_types.h"
#include "BKE_object.h"
#include "BKE_displist.h"
#include "BKE_ama.h"
#include "BKE_anim.h"
#include "BKE_cdderivedmesh.h"
#include "BKE_mesh.h"
#include "BKE_object.h"

/* XXX This function fixes bad merging code, in some cases removing vertices creates indices > maxvert */
int test_index_face_maxvert(MFace *mface, CustomData *fdata, int mfindex, int nr, int maxvert)
{
	if(mface->v1 >= maxvert) {
		// printf("bad index in array\n");
		mface->v1= maxvert - 1;
	}
	if(mface->v2 >= maxvert) {
		// printf("bad index in array\n");
		mface->v2= maxvert - 1;
	}
	if(mface->v3 >= maxvert) {
		// printf("bad index in array\n");
		mface->v3= maxvert - 1;
	}
	if(mface->v4 >= maxvert) {
		// printf("bad index in array\n");
		mface->v4= maxvert - 1;
	}

	return test_index_face(mface, fdata, mfindex, nr);
}


/* indexMap - an array of IndexMap entries
 * oldIndex - the old index to map
 * copyNum - the copy number to map to (original = 0, first copy = 1, etc.)
 */
int calc_mapping(IndexMapEntry *indexMap, int oldIndex, int copyNum)
{
	if(indexMap[oldIndex].merge < 0) {
		/* vert wasn't merged, so use copy of this vert */
		return indexMap[oldIndex].new + copyNum;
	} else if(indexMap[oldIndex].merge == oldIndex) {
		/* vert was merged with itself */
		return indexMap[oldIndex].new;
	} else {
		/* vert was merged with another vert */
		/* follow the chain of merges to the end, or until we've passed
		* a number of vertices equal to the copy number
		*/
		if(copyNum <= 0)
			return indexMap[oldIndex].new;
		else
			return calc_mapping(indexMap, indexMap[oldIndex].merge,
						copyNum - 1);
	}
}


float length_fitcurve(ArrayModifierData *amd, struct Scene *scene)
{
	float length = 0;

	Curve *cu = amd->curve_ob->data;
	if(cu) {
		float tmp_mat[3][3];
		float scale;

		BKE_object_to_mat3(amd->curve_ob, tmp_mat);
		scale = mat3_to_scale(tmp_mat);

		if(!cu->path) {
			cu->flag |= CU_PATH; // needed for path & bevlist
			BKE_displist_make_curveTypes(scene, amd->curve_ob, 0);
		}
		if(cu->path)
			length = scale*cu->path->totdist;
	}
	return length;
}


int length_to_count(float length, const float offset[3])
{
	int count = 0;

	float dist = sqrt(dot_v3v3(offset, offset));

	if(dist > 1e-6f)
		/* this gives length = first copy start to last copy end
		add a tiny offset for floating point rounding errors */
		count = (length + 1e-6f) / dist;
	else
		/* if the offset has no translation, just make one copy */
		count = 1;
	return count;
}


float count_to_length(int count, const float offset[3])
{
	float length = 0;
	float dist = sqrt(dot_v3v3(offset, offset));

	if(dist > 1e-6f)
		/* this gives length = first copy start to last copy end
		add a tiny offset for floating point rounding errors */
		//count = (length + 1e-6f) / dist;
		length = count * dist - 1e-6f;
	else
		/* if the offset has no translation, just make one copy */
		length = 1;
	return length;
}


//generates a psuedo-random float between 0.0 and max
float f_rand_max(float max)
{
	return BLI_frand()*max;
}


void array_offset(const float max_off[3], float rit[3], int prop, int sign, int seed)
{
	int j;

	BLI_srandom(seed);

	rit[0] = f_rand_max(max_off[0]);
	if (sign & MOD_ARR_SIGN_L) {
		if (sign & MOD_ARR_SIGN_P) {
			j = BLI_rand() % 2;
			if (j == 0)
				rit[0] = rit[0]*(-1);
		}
		else
			rit[0] = rit[0]*(-1);
	}

	if (!(prop & MOD_ARR_PROP)) {
		rit[1] = f_rand_max(max_off[1]);
		if (sign & MOD_ARR_SIGN_L) {
			if (sign & MOD_ARR_SIGN_P) {
				j = BLI_rand() % 2;
				if (j == 0)
					rit[1] = rit[1]*(-1);
			}
			else
				rit[1] = rit[1]*(-1);
		}

		rit[2] = f_rand_max(max_off[2]);
		if (sign & MOD_ARR_SIGN_L) {
			if (sign & MOD_ARR_SIGN_P) {
				j = BLI_rand() % 2;
				if (j == 0)
					rit[2] = rit[2]*(-1);
			}
			else
				rit[2] = rit[2]*(-1);
		}
	}
	else {
		rit[1] = rit[0];
		rit[2] = rit[0];
	}
}

void init_mat_oc(const int start, const int end, int *vet_mc)
{
	int i;

	for (i=start; i<end; i++)
		vet_mc[i] = 0;
}


void init_offset(const int start, const int end, ArrayModifierData *ar)
{
	int i;

	for (i=start; i<end; i++) {
		zero_v3(ar->Mem_Ob[i].rot);
		ar->Mem_Ob[i].scale[0] = ar->Mem_Ob[i].scale[1] = ar->Mem_Ob[i].scale[2] = 1;
		zero_v3(ar->Mem_Ob[i].loc);
		zero_v3(ar->Mem_Ob[i].cu_cent);
		zero_v4(ar->Mem_Ob[i].cu_loc);
		ar->Mem_Ob[i].id_mat = 0;
		ar->Mem_Ob[i].transform = 0;
		ar->Mem_Ob[i].rand_group_obj = 0;
	}
}


void create_offset(const int n, const int totmat, ArrayModifierData *ar, Object *ob)
{
	float loc[3];
	float rot[3];
	float scale[3];
	int i, seed_t, seed_g, seed_m, act_mat = 0;
	int cont_mat = ar->cont_mat-1;
	Group *group = NULL;

	if(ob->dup_group!=NULL)
		group= ob->dup_group;

	seed_t = ar->seed[0];
	seed_g = ar->seed[1];
	seed_m = ar->seed[2];
	for (i=0; i < n-1; i++) {

		loc[0]=loc[1]=loc[2]=0;
		rot[0]=rot[1]=rot[2]=0;
		scale[0]=scale[1]=scale[2]=1;

		if (ar->mode & MOD_ARR_MOD_ADV) {
			if (ar->lock & MOD_ARR_LOCK_ROT) {
				if ((ar->rot_offset[0]!=0) || (ar->rot_offset[1]!=0) || (ar->rot_offset[2]!=0)) {
					array_offset(ar->rot_offset, rot, !MOD_ARR_PROP, ar->sign, seed_t);
					ar->Mem_Ob[i].transform = 1;
				}
			}
			if (ar->lock & MOD_ARR_LOCK_SCALE) {
				if ((ar->scale_offset[0]!=1) || (ar->scale_offset[1]!=1) || (ar->scale_offset[2]!=1)) {
					array_offset(ar->scale_offset, scale, ar->flag_offset, ar->sign, seed_t);
					ar->Mem_Ob[i].transform = 1;
				}
			}
			if (ar->lock & MOD_ARR_LOCK_LOC) {
				if ((ar->loc_offset[0]!=0) || (ar->loc_offset[1]!=0) || (ar->loc_offset[2]!=0)) {
					array_offset(ar->loc_offset, loc, !MOD_ARR_PROP, ar->sign, seed_t);
					ar->Mem_Ob[i].transform = 1;
				}
			}
			if (ar->Mem_Ob[i].transform) {
				/* Rotation */
				if (ar->lock & MOD_ARR_LOCK_ROT)
					copy_v3_v3(ar->Mem_Ob[i].rot, rot);
				/* Scaling */
				if (ar->lock & MOD_ARR_LOCK_SCALE)
					copy_v3_v3(ar->Mem_Ob[i].scale, scale);
				/* Location */
				if (ar->lock & MOD_ARR_LOCK_LOC)
					copy_v3_v3(ar->Mem_Ob[i].loc, loc);
			}
		}
		if (ar->mode & MOD_ARR_MOD_ADV_CLONE) {
			if(group!=NULL) {
				if (ar->rand_group & MOD_ARR_RAND_GROUP) {
					BLI_srandom(seed_g);
					ar->Mem_Ob[i].rand_group_obj = BLI_rand() % BLI_countlist(&group->gobject);
					ar->Mem_Ob[i].rand_group_obj++;
				}
			}
		}
		/* Random Material Clone */
		if (ar->mode & MOD_ARR_MOD_ADV_MAT) {
			if (totmat>1) {
				/* Random */
				if ((ar->rand_mat & MOD_ARR_MAT) && (ar->mat_ob & MOD_ARR_AR_MAT_RND)) {
					BLI_srandom(seed_m);
					ar->Mem_Ob[i].id_mat = BLI_rand() % totmat;
				}
				else { /* Sequence */
					if (cont_mat == 0 ){
						cont_mat = ar->cont_mat;
						if (act_mat + 1 < totmat)
							act_mat++;
						else
							act_mat = 0;
					}
					ar->Mem_Ob[i].id_mat = act_mat;
					cont_mat--;
				}
			}
		}
		seed_t++;
		seed_g++;
		seed_m++;
	}

	/* Random Material Object Cap */
	if (ar->mode & MOD_ARR_MOD_ADV_MAT) {
		if (totmat>1) {
			if ((ar->rand_mat & MOD_ARR_MAT) && (ar->mat_ob & MOD_ARR_SC_MAT_RND)) {
				ar->Mem_Mat_Ob.start_cap = BLI_rand() % totmat;
			}
			if ((ar->rand_mat & MOD_ARR_MAT) && (ar->mat_ob & MOD_ARR_MC_MAT_RND)) {
				int dim_mc;
				if (ar->mode & MOD_ARR_MOD_ADV_MID)
					dim_mc = ar->count_mc;
				else
					dim_mc = 1;
				for (i=0; i < dim_mc; i++) {
					ar->Mem_Mat_Ob.mid_cap[i] = BLI_rand() % totmat;
				}
			}
			if ((ar->rand_mat & MOD_ARR_MAT) && (ar->mat_ob & MOD_ARR_EC_MAT_RND)) {
				ar->Mem_Mat_Ob.end_cap = BLI_rand() % totmat;
			}
		}
	}
}


static void init_curve_deform(Object *par, Object *ob, CurveDeform *cd, int dloc)
{
	invert_m4_m4(ob->imat, ob->obmat);
	mul_m4_m4m4(cd->objectspace, ob->imat, par->obmat);
	invert_m4_m4(cd->curvespace, cd->objectspace);
	copy_m3_m4(cd->objectspace3, cd->objectspace);

	// offset vector for 'no smear'
	if(dloc) {
		invert_m4_m4(par->imat, par->obmat);
		mul_v3_m4v3(cd->dloc, par->imat, ob->obmat[3]);
	}
	else {
		cd->dloc[0]=cd->dloc[1]=cd->dloc[2]= 0.0f;
	}

	cd->no_rot_axis= 0;
}


/* this makes sure we can extend for non-cyclic. *vec needs 4 items! */
static int where_on_path_deform(Object *ob, float ctime, float *vec, float *dir, float *quat, float *radius)	/* returns OK */
{
	Curve *cu= ob->data;
	BevList *bl;
	float ctime1;
	int cycl=0;

	/* test for cyclic */
	bl= cu->bev.first;
	if (!bl->nr) return 0;
	if(bl && bl->poly> -1) cycl= 1;

	if(cycl==0) {
		ctime1= CLAMPIS(ctime, 0.0f, 1.0f);
	}
	else ctime1= ctime;

	/* vec needs 4 items */
	if(where_on_path(ob, ctime1, vec, dir, quat, radius, NULL)) {

		if(cycl==0) {
			Path *path= cu->path;
			float dvec[3];

			if(ctime < 0.0f) {
				sub_v3_v3v3(dvec, path->data[1].vec, path->data[0].vec);
				mul_v3_fl(dvec, ctime*(float)path->len);
				add_v3_v3(vec, dvec);
				if(quat) copy_qt_qt(quat, path->data[0].quat);
				if(radius) *radius= path->data[0].radius;
				print_v4("Vec", vec);
			}
			else if(ctime > 1.0f) {
				sub_v3_v3v3(dvec, path->data[path->len-1].vec, path->data[path->len-2].vec);
				mul_v3_fl(dvec, (ctime-1.0f)*(float)path->len);
				add_v3_v3(vec, dvec);
				if(quat) copy_qt_qt(quat, path->data[path->len-1].quat);
				if(radius) *radius= path->data[path->len-1].radius;
				/* weight - not used but could be added */
			}

		}
		return 1;
	}
	return 0;
}


	/* for each point, rotate & translate to curve */
	/* use path, since it has constant distances */
	/* co: local coord, result local too */
	/* returns quaternion for rotation, using cd->no_rot_axis */
	/* axis is using another define!!! */
static int calc_curve_deform(Scene *scene, Object *par, float *co, CurveDeform *cd, float *loc, float *cent)
{
	Curve *cu= par->data;
	float fac, dir[3], new_quat[4], radius;
	short /*upflag, */ index, axis = 1;

	index= axis-1;
	if(index>2)
		index -= 3; /* negative  */

	/* to be sure, mostly after file load */
	if(cu->path==NULL) {
		BKE_displist_make_curveTypes(scene, par, 0);
		if(cu->path==NULL) return 0;	// happens on append...
	}

	if(ELEM3(axis, OB_NEGX+1, OB_NEGY+1, OB_NEGZ+1)) { /* OB_NEG# 0-5, MOD_CURVE_POS# 1-6 */
		if(cu->flag & CU_STRETCH)
			fac= (-co[index]-cd->dmax[index])/(cd->dmax[index] - cd->dmin[index]);
		else
			fac= (cd->dloc[index])/(cu->path->totdist) - (co[index]-cd->dmax[index])/(cu->path->totdist);
	}
	else {
		if(cu->flag & CU_STRETCH)
			fac= (co[index]-cd->dmin[index])/(cd->dmax[index] - cd->dmin[index]);
		else
			fac= (cd->dloc[index])/(cu->path->totdist) + (co[index]-cd->dmin[index])/(cu->path->totdist);
	}

	if( where_on_path_deform(par, fac, loc, dir, new_quat, &radius)) {	/* returns OK */
		float quat[4]; //, cent[3];
		if(cd->no_rot_axis) {	/* set by caller */

			/* this is not exactly the same as 2.4x, since the axis is having rotation removed rather than
			 * changing the axis before calculating the tilt but serves much the same purpose */
			float dir_flat[3]={0,0,0}, q[4];
			copy_v3_v3(dir_flat, dir);
			dir_flat[cd->no_rot_axis-1]= 0.0f;

			normalize_v3(dir);
			normalize_v3(dir_flat);

			rotation_between_vecs_to_quat(q, dir, dir_flat); /* Could this be done faster? */

			mul_qt_qtqt(new_quat, q, new_quat);
		}


		/* Logic for 'cent' orientation *
		 *
		 * The way 'co' is copied to 'cent' may seem to have no meaning, but it does.
		 *
		 * Use a curve modifier to stretch a cube out, color each side RGB, positive side light, negative dark.
		 * view with X up (default), from the angle that you can see 3 faces RGB colors (light), anti-clockwise
		 * Notice X,Y,Z Up all have light colors and each ordered CCW.
		 *
		 * Now for Neg Up XYZ, the colors are all dark, and ordered clockwise - Campbell
		 *
		 * note: moved functions into quat_apply_track/vec_apply_track
		 * */
		copy_qt_qt(quat, new_quat);
		copy_v3_v3(cent, co);

		/* zero the axis which is not used,
		 * the big block of text above now applies to these 3 lines */
		quat_apply_track(quat, axis-1, (axis==1 || axis==3) ? 1:0); /* up flag is a dummy, set so no rotation is done */
		vec_apply_track(cent, axis-1);
		cent[axis < 4 ? axis-1 : axis-4]= 0.0f;


		/* scale if enabled */
		if(cu->flag & CU_PATH_RADIUS)
			mul_v3_fl(cent, radius);

		/* local rotation */
		normalize_qt(quat);
		mul_qt_v3(quat, cent);

		/* translation */
		/*add_v3_v3v3(co, cent, loc);*/

		return 1;
	}
	return 0;
}


void array_to_curve(Scene *scene, Object *cuOb, Object *target, float *vertexCos, float *vec, float *quat)
{
	Curve *cu;
	CurveDeform cd;

	if(cuOb->type != OB_CURVE)
		return;
	cu = cuOb->data;
	cu->flag |= (CU_PATH|CU_FOLLOW);
	init_curve_deform(cuOb, target, &cd, (cu->flag & CU_STRETCH)==0);
	calc_curve_deform(scene, cuOb, vertexCos, &cd, vec, quat);
}


/** Divide the line segments associated with the currently selected
 * curve nodes (Bezier or NURB). If there are no valid segment
 * selections within the current selection, nothing happens.
 */
void dist_point_nurb(Object *obedit, int number_cuts)
{
	Curve *cu= obedit->data;
	EditNurb *editnurb= cu->editnurb;
	Nurb *nu;
	BezTriple *prevbezt, *bezt, *beztnew, *beztn;
	BPoint *bp, *prevbp, *bpnew, *bpn;
	float vec[15];
	int a, amount, i;
	float factor;

	for(nu= editnurb->nurbs.first; nu; nu= nu->next) {
		amount= 0;
		if(nu->type == CU_BEZIER) {
		/*
		   Insert a point into a 2D Bezier curve.
		   Endpoints are preserved. Otherwise, all selected and inserted points are
		   newly created. Old points are discarded.
		*/
			/* count */
			if(nu->flagu & CU_NURB_CYCLIC) {
				a = nu->pntsu;
				bezt = nu->bezt;
				prevbezt = bezt+(a-1);
			}
			else {
				a = nu->pntsu-1;
				prevbezt = nu->bezt;
				bezt = prevbezt+1;
			}
			/*while(a--) {
				if( BEZSELECTED_HIDDENHANDLES(cu, prevbezt) && BEZSELECTED_HIDDENHANDLES(cu, bezt) ) amount+=number_cuts;
				prevbezt= bezt;
				bezt++;
			}*/

			if(amount) {
				/* insert */
				beztnew =
					(BezTriple*)MEM_mallocN((amount + nu->pntsu) * sizeof(BezTriple), "subdivNurb");
				/*beztn= beztnew;
				if(nu->flagu & CU_NURB_CYCLIC) {
					a = nu->pntsu;
					bezt = nu->bezt;
					prevbezt = bezt+(a-1);
				}
				else {
					a = nu->pntsu-1;
					prevbezt = nu->bezt;
					bezt = prevbezt+1;
				}*/
				while(a--) {
					memcpy(beztn, prevbezt, sizeof(BezTriple));
					//keyIndex_updateBezt(editnurb, prevbezt, beztn, 1);
					beztn++;

					if( BEZSELECTED_HIDDENHANDLES(cu, prevbezt) && BEZSELECTED_HIDDENHANDLES(cu, bezt) ) {
						float prevvec[3][3];

						memcpy(prevvec, prevbezt->vec, sizeof(float) * 9);

						for (i = 0; i < number_cuts; i++) {
							factor = 1.0f / (number_cuts + 1 - i);

							memcpy(beztn, bezt, sizeof(BezTriple));

							/* midpoint subdividing */
							interp_v3_v3v3(vec, prevvec[1], prevvec[2], factor);
							interp_v3_v3v3(vec+3, prevvec[2], bezt->vec[0], factor);
							interp_v3_v3v3(vec+6, bezt->vec[0], bezt->vec[1], factor);

							interp_v3_v3v3(vec+9, vec, vec+3, factor);
							interp_v3_v3v3(vec+12, vec+3, vec+6, factor);

							/* change handle of prev beztn */
							copy_v3_v3((beztn-1)->vec[2], vec);
							/* new point */
							copy_v3_v3(beztn->vec[0], vec+9);
							interp_v3_v3v3(beztn->vec[1], vec+9, vec+12, factor);
							copy_v3_v3(beztn->vec[2], vec+12);
							/* handle of next bezt */
							if(a==0 && i == number_cuts - 1 && (nu->flagu & CU_NURB_CYCLIC)) {copy_v3_v3(beztnew->vec[0], vec+6);}
							else {copy_v3_v3(bezt->vec[0], vec+6);}

							beztn->radius = (prevbezt->radius + bezt->radius)/2;
							beztn->weight = (prevbezt->weight + bezt->weight)/2;

							memcpy(prevvec, beztn->vec, sizeof(float) * 9);
							beztn++;
						}
					}

					prevbezt= bezt;
					bezt++;
				}
				/* last point */
				if((nu->flagu & CU_NURB_CYCLIC)==0) {
					memcpy(beztn, prevbezt, sizeof(BezTriple));
					//keyIndex_updateBezt(editnurb, prevbezt, beztn, 1);
				}

				MEM_freeN(nu->bezt);
				nu->bezt= beztnew;
				nu->pntsu+= amount;

				//calchandlesNurb(nu);
			}
		} /* End of 'if(nu->type == CU_BEZIER)' */
		else if (nu->pntsv==1) {
		/*
		   All flat lines (ie. co-planar), except flat Nurbs. Flat NURB curves
		   are handled together with the regular NURB plane division, as it
		   should be. I split it off just now, let's see if it is
		   stable... nzc 30-5-'00
		 */
			/* count */
			if(nu->flagu & CU_NURB_CYCLIC) {
				a= nu->pntsu;
				bp= nu->bp;
				prevbp= bp+(a-1);
			}
			else {
				a= nu->pntsu-1;
				prevbp= nu->bp;
				bp= prevbp+1;
			}
			while(a--) {
				if( (bp->f1 & SELECT) && (prevbp->f1 & SELECT) ) amount+=number_cuts;
				prevbp= bp;
				bp++;
			}

			if(amount) {
				/* insert */
				bpnew =
					(BPoint*)MEM_mallocN((amount + nu->pntsu) * sizeof(BPoint), "subdivNurb2");
				bpn= bpnew;

				if(nu->flagu & CU_NURB_CYCLIC) {
					a= nu->pntsu;
					bp= nu->bp;
					prevbp= bp+(a-1);
				}
				else {
					a= nu->pntsu-1;
					prevbp= nu->bp;
					bp= prevbp+1;
				}
				while(a--) {
					memcpy(bpn, prevbp, sizeof(BPoint));
					//keyIndex_updateBP(editnurb, prevbp, bpn, 1);
					bpn++;

					if( (bp->f1 & SELECT) && (prevbp->f1 & SELECT) ) {
				 // printf("*** subdivideNurb: insert 'linear' point\n");
						for (i = 0; i < number_cuts; i++) {
							factor = (float)(i + 1) / (number_cuts + 1);

							memcpy(bpn, bp, sizeof(BPoint));
							interp_v4_v4v4(bpn->vec, prevbp->vec, bp->vec, factor);
							bpn++;
						}

					}
					prevbp= bp;
					bp++;
				}
				if((nu->flagu & CU_NURB_CYCLIC)==0) { /* last point */
					memcpy(bpn, prevbp, sizeof(BPoint));
					//keyIndex_updateBP(editnurb, prevbp, bpn, 1);
				}

				MEM_freeN(nu->bp);
				nu->bp= bpnew;
				nu->pntsu+= amount;

				if(nu->type & CU_NURBS) {
					//nurbs_knot_calc_u(nu);
				}
			}
		} /* End of 'else if(nu->pntsv==1)' */
	}
}

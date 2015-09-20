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
 * The Original Code is Copyright (C) 2015 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Blender Foundation,
 *                 Bastien Montagne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __EIGEN3_EIGENVALUES_C_API_H__
#define __EIGEN3_EIGENVALUES_C_API_H__

#ifdef __cplusplus
extern "C" {
#endif

bool EG3_self_adjoint_eigen_solve(const int size, const float *matrix, float *r_eigen_values, float *r_eigen_vectors);

#ifdef __cplusplus
}
#endif

#endif  /* __EIGEN3_EIGENVALUES_C_API_H__ */

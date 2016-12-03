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
 * Contributor(s): Esteban Tovagliari, Cedric Paille, Kevin Dietrich
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include "abc_transform.h"

#include <OpenEXR/ImathBoxAlgo.h>

#include "abc_util.h"

extern "C" {
#include "DNA_object_types.h"

#include "BLI_math.h"

#include "BKE_object.h"
}

using Alembic::AbcGeom::OObject;
using Alembic::AbcGeom::OXform;

/* ************************************************************************** */

static bool has_parent_camera(Object *ob)
{
	if (!ob->parent) {
		return false;
	}

	Object *parent = ob->parent;

	if (parent->type == OB_CAMERA) {
		return true;
	}

	return has_parent_camera(parent);
}

/* ************************************************************************** */

AbcTransformWriter::AbcTransformWriter(Object *ob,
                                       const OObject &abc_parent,
                                       AbcTransformWriter *parent,
                                       unsigned int time_sampling,
                                       ExportSettings &settings)
    : AbcObjectWriter(NULL, ob, time_sampling, settings, parent)
{
	m_is_animated = hasAnimation(m_object);
	m_parent = NULL;

	if (!m_is_animated) {
		time_sampling = 0;
	}

	m_xform = OXform(abc_parent, get_id_name(m_object), time_sampling);
	m_schema = m_xform.getSchema();
}

void AbcTransformWriter::do_write()
{
	if (m_first_frame) {
		m_visibility = Alembic::AbcGeom::CreateVisibilityProperty(m_xform, m_xform.getSchema().getTimeSampling());
	}

	m_visibility.set(!(m_object->restrictflag & OB_RESTRICT_VIEW));

	if (!m_first_frame && !m_is_animated) {
		return;
	}

	float mat[4][4];
	create_transform_matrix(m_object, mat);

	/* Only apply rotation to root camera, parenting will propagate it. */
	if (m_object->type == OB_CAMERA && !has_parent_camera(m_object)) {
		float rot_mat[4][4];
		axis_angle_to_mat4_single(rot_mat, 'X', -M_PI_2);
		mul_m4_m4m4(mat, mat, rot_mat);
	}

	if (!m_object->parent) {
		/* Only apply scaling to root objects, parenting will propagate it. */
		float scale_mat[4][4];
		scale_m4_fl(scale_mat, m_settings.global_scale);
		mul_m4_m4m4(mat, mat, scale_mat);
		mul_v3_fl(mat[3], m_settings.global_scale);
	}

	m_matrix = convert_matrix(mat);

	m_sample.setMatrix(m_matrix);
	m_schema.set(m_sample);
}

Imath::Box3d AbcTransformWriter::bounds()
{
	Imath::Box3d bounds;

	for (int i = 0; i < m_children.size(); ++i) {
		Imath::Box3d box(m_children[i]->bounds());
		bounds.extendBy(box);
	}

	return Imath::transform(bounds, m_matrix);
}

bool AbcTransformWriter::hasAnimation(Object */*ob*/) const
{
	/* TODO(kevin): implement this. */
	return true;
}

/* ************************************************************************** */

AbcEmptyReader::AbcEmptyReader(const Alembic::Abc::IObject &object, ImportSettings &settings)
    : AbcObjectReader(object, settings)
{
	Alembic::AbcGeom::IXform xform(object, Alembic::AbcGeom::kWrapExisting);
	m_schema = xform.getSchema();

	get_min_max_time(m_iobject, m_schema, m_min_time, m_max_time);
}

bool AbcEmptyReader::valid() const
{
	return m_schema.valid();
}

void AbcEmptyReader::readObjectData(Main *bmain, float /*time*/)
{
	m_object = BKE_object_add_only_object(bmain, OB_EMPTY, m_object_name.c_str());
	m_object->data = NULL;
}

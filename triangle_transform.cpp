#include "triangle_transform.h"

#include "core/error/error_macros.h"
#include "core/math/basis.h"
#include "core/math/transform_3d.h"
#include "core/math/vector3.h"
#include "core/object/class_db.h"
#include "core/typedefs.h"
#include "core/variant/variant.h"
#include "pseudo_inverse.h"

TriangleTransform::TriangleTransform(const PackedVector3Array &vertices, const Transform3D &point_tf)
{
	// Compute vertice lengths
	if(unlikely(vertices.size() != 3))
	{
		this->coord_lengths = Vector3(0, 0, 0);
	}
	else
	{
		this->update_transform(vertices[0], vertices[1], vertices[2], point_tf);
	}
}

TriangleTransform::TriangleTransform(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
                                     const Transform3D &point_tf)
{
	this->update_transform(vert0, vert1, vert2, point_tf);
}

Vector3 TriangleTransform::adjust_point(const PackedVector3Array &vertices, const float &normal_offset) const
{
	ERR_FAIL_COND_V(vertices.size() != 3, Vector3(0, 0, 0));
	return this->adjust_point(vertices[0], vertices[1], vertices[2], normal_offset);
}

Vector3 TriangleTransform::adjust_point(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
                                        const float &normal_offset) const
{
	// Create coordinate system
	Vector3 dir1       = vert1 - vert0;
	Vector3 dir2       = vert2 - vert0;
	const Vector3 norm = compute_norm(dir1, dir2);

	dir1 /= dir1.dot(dir1);
	dir2 /= dir2.dot(dir2);

	return vert0 + this->compute_rel_point_from_system(dir1, dir2, norm, normal_offset);
}

Transform3D TriangleTransform::adjust_transform(const PackedVector3Array &vertices, const float &normal_offset)
{
	ERR_FAIL_COND_V(vertices.size() != 3, Transform3D());
	return this->adjust_transform(vertices[0], vertices[1], vertices[2], normal_offset);
}

Transform3D TriangleTransform::adjust_transform(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
                                                const float &normal_offset)
{
	Transform3D ret;

	// Create coordinate system
	Vector3 dir1       = vert1 - vert0;
	Vector3 dir2       = vert2 - vert0;
	const Vector3 norm = compute_norm(dir1, dir2);

	dir1 /= dir1.dot(dir1);
	dir2 /= dir2.dot(dir2);

	// Compute updated point
	ret.origin = vert0 + this->compute_rel_point_from_system(dir1, dir2, norm, normal_offset);

	// Compute updated rotation
	ret.basis = this->compute_normal_rotation(norm, dir1);

	return ret;
}

void TriangleTransform::update_transform(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
                                         const Transform3D &point_tf)
{
	// Create coordinate system
	const Vector3 dir1 = vert1 - vert0;
	const Vector3 dir2 = vert2 - vert0;
	const Vector3 norm = compute_norm(dir1, dir2);

	this->coord_lengths     = compute_coord_lengths_from_system(dir1, dir2, norm, point_tf.origin - vert0);
	this->original_rotation = Basis::looking_at(norm, dir1);
}

Vector3 TriangleTransform::compute_coord_lengths(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
                                                 const Transform3D &point_tf)
{
	// Create coordinate system
	const Vector3 dir1 = vert1 - vert0;
	const Vector3 dir2 = vert2 - vert0;
	const Vector3 norm = compute_norm(dir1, dir2);

	// Compute vertice lengths
	const Vector3 rel_point = point_tf.origin - vert0;
	return compute_coord_lengths_from_system(dir1, dir2, norm, rel_point);
}

void TriangleTransform::_bind_methods()
{
	{
		using adjust_point_fcn_t = Vector3 (TriangleTransform::*)(const PackedVector3Array &, const float &) const;
		ClassDB::bind_method(D_METHOD("adjust_point", "vertices", "normal_offset"),
		                     (adjust_point_fcn_t)&TriangleTransform::adjust_point);
	}

	{
		using adjust_point_fcn_t =
			Vector3 (TriangleTransform::*)(const Vector3 &, const Vector3 &, const Vector3 &, const float &) const;
		ClassDB::bind_method(D_METHOD("adjust_point_from_verts", "vert0", "vert1", "vert2", "normal_offset"),
		                     (adjust_point_fcn_t)&TriangleTransform::adjust_point);
	}

	{
		using adjust_transform_fcn_t = Transform3D (TriangleTransform::*)(const PackedVector3Array &, const float &);
		ClassDB::bind_method(D_METHOD("adjust_transform", "vertices", "normal_offset"),
		                     (adjust_transform_fcn_t)&TriangleTransform::adjust_transform);
	}

	{
		using adjust_transform_fcn_t =
			Transform3D (TriangleTransform::*)(const Vector3 &, const Vector3 &, const Vector3 &, const float &);
		ClassDB::bind_method(D_METHOD("adjust_transform_from_verts", "vert0", "vert1", "vert2", "normal_offset"),
		                     (adjust_transform_fcn_t)&TriangleTransform::adjust_transform);
	}

	ClassDB::bind_method(D_METHOD("get_coord_lengths"), &TriangleTransform::get_coord_lengths);
	ClassDB::bind_method(D_METHOD("set_coord_lengths", "coord_lengths"), &TriangleTransform::set_coord_lengths);
	ClassDB::add_property("TriangleTransform", PropertyInfo(Variant::Type::VECTOR3, "coord_lengths"),
	                      "set_coord_lengths", "get_coord_lengths");

	ClassDB::bind_method(D_METHOD("get_original_rotation"), &TriangleTransform::get_original_rotation);
}

Vector3 TriangleTransform::compute_coord_lengths_from_system(const Vector3 &dir1, const Vector3 &dir2,
                                                             const Vector3 &norm, const Vector3 &rel_point)
{
	// Create coordinate system
	PackedVector3Array coord_matrix;
	coord_matrix.resize(3);
	*(coord_matrix.ptrw() + 0) = dir1 / dir1.dot(dir1);
	*(coord_matrix.ptrw() + 1) = dir2 / dir2.dot(dir2);
	*(coord_matrix.ptrw() + 2) = norm;

	// Compute vertice lengths
	return pseudo_inverse_mult(coord_matrix, rel_point);
}

Vector3 TriangleTransform::compute_rel_point_from_system(const Vector3 &dir1_adj, const Vector3 &dir2_adj,
                                                         const Vector3 &norm, const float &normal_offset) const
{
	return Vector3(dir1_adj[0] * this->coord_lengths[0] + dir2_adj[0] * this->coord_lengths[1] +
	                   norm[0] * this->coord_lengths[2],
	               dir1_adj[1] * this->coord_lengths[0] + dir2_adj[1] * this->coord_lengths[1] +
	                   norm[1] * this->coord_lengths[2],
	               dir1_adj[2] * this->coord_lengths[0] + dir2_adj[2] * this->coord_lengths[1] +
	                   norm[2] * this->coord_lengths[2]) +
	       normal_offset * norm;
}

Basis TriangleTransform::compute_normal_rotation(const Vector3 &normal, const Vector3 &dir1)
{
	return Basis::looking_at(normal, dir1); //*this->original_rotation.inverse();
}

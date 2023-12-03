#pragma once

#include "core/math/transform_3d.h"
#include "core/math/vector3.h"
#include "core/object/object.h"
#include "core/variant/variant.h"
#include <core/object/ref_counted.h>

class TriangleTransform : public RefCounted
{
	GDCLASS(TriangleTransform, RefCounted);

	public:
	TriangleTransform() = default;
	TriangleTransform(const PackedVector3Array &vertices, const Transform3D &point_tf);
	TriangleTransform(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2, const Transform3D &point_tf);

	Vector3 adjust_point(const PackedVector3Array &vertices, const float &normal_offset = 0.0f) const;
	Vector3 adjust_point(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
	                     const float &normal_offset = 0.0f) const;

	Transform3D adjust_transform(const PackedVector3Array &vertices, const float &normal_offset = 0.0f);
	Transform3D adjust_transform(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
	                             const float &normal_offset = 0.0f);

	void update_transform(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
	                      const Transform3D &point_tf);

	static Vector3 compute_coord_lengths(const Vector3 &vert0, const Vector3 &vert1, const Vector3 &vert2,
	                                     const Transform3D &point_tf);

	protected:
	static void _bind_methods();

	private:
	// Length of coordinate system parts:
	// coord_lengths[0]: Length along vertices[1]-vertices[0]
	// coord_lengths[1]: Length along vertices[2]-vertices[0]
	// coord_lengths[2]: Length along vertice normal
	Vector3 coord_lengths;

	// Original rotation on construction
	Basis original_rotation;

	Vector3 get_coord_lengths()
	{
		return this->coord_lengths;
	}

	void set_coord_lengths(const Vector3 &_coord_lengths)
	{
		this->coord_lengths = _coord_lengths;
	}

	static Vector3 compute_coord_lengths_from_system(const Vector3 &dir1, const Vector3 &dir2, const Vector3 &norm,
	                                                 const Vector3 &rel_point);

	Vector3 compute_rel_point_from_system(const Vector3 &dir1_adj, const Vector3 &dir2_adj, const Vector3 &norm,
	                                      const float &normal_offset = 0.0f) const;
	Basis compute_normal_rotation(const Vector3 &normal, const Vector3 &dir1);

	static inline Vector3 compute_norm(const Vector3 &dir1, const Vector3 &dir2)
	{
		return dir2.cross(dir1).normalized();
	}
};

#pragma once

#include <core/object/ref_counted.h>
#include <scene/3d/mesh_instance_3d.h>

struct MeshTrianglePoint : public RefCounted
{
	static constexpr auto INVALID_ID   = std::numeric_limits<uint32_t>::max();
	static constexpr auto INVALID_DIST = std::numeric_limits<float>::infinity();

	MeshInstance3D *mesh_instance = nullptr;
	uint32_t surface_id           = INVALID_ID;
	uint32_t index_id             = INVALID_ID;
	float ray_origin_dist         = INVALID_DIST;

	int32_t vertex_ids[3];

	Vector3 point_on_triangle;

	GDCLASS(MeshTrianglePoint, RefCounted);

	public:
	MeshTrianglePoint() = default;
	MeshTrianglePoint(MeshInstance3D *_mesh_instance, uint32_t _surface_id, uint32_t _index_id, float _ray_origin_dist);

	protected:
	static void _bind_methods();

	MeshInstance3D *get_mesh_instance()
	{
		return this->mesh_instance;
	}

	void set_mesh_instance(MeshInstance3D *_mesh_instance)
	{
		this->mesh_instance = _mesh_instance;
	}

	uint32_t get_surface_id()
	{
		return this->surface_id;
	}

	void set_surface_id(uint32_t _surface_id)
	{
		this->surface_id = _surface_id;
	}

	uint32_t get_index_id()
	{
		return this->index_id;
	}

	void set_index_id(uint32_t _index_id)
	{
		this->index_id = _index_id;
	}

	float get_ray_origin_dist()
	{
		return this->ray_origin_dist;
	}

	void set_ray_origin_dist(float _ray_origin_dist)
	{
		this->ray_origin_dist = _ray_origin_dist;
	}

	PackedInt32Array get_vertex_ids()
	{
		return PackedInt32Array({this->vertex_ids[0], this->vertex_ids[1], this->vertex_ids[2]});
	}

	void set_vertex_ids(const PackedInt32Array &_vertex_ids)
	{
		if(_vertex_ids.size() != 3)
		{
			return;
		}

		this->vertex_ids[0] = _vertex_ids[0];
		this->vertex_ids[1] = _vertex_ids[1];
		this->vertex_ids[2] = _vertex_ids[2];
	}

	Vector3 get_point_on_triangle()
	{
		return this->point_on_triangle;
	}

	void set_point_on_triangle(Vector3 _point_on_triangle)
	{
		this->point_on_triangle = _point_on_triangle;
	}
};

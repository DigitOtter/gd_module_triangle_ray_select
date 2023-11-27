#pragma once

#include "core/math/math_defs.h"
#include "core/variant/variant.h"
#include "shaders/triangle_ray_select.glsl.gen.h"

#include <core/object/ref_counted.h>
#include <scene/3d/camera_3d.h>
#include <scene/3d/mesh_instance_3d.h>
#include <servers/rendering/renderer_rd/storage_rd/mesh_storage.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>

struct MeshSurfaceIndex : public RefCounted
{
	static constexpr auto INVALID_ID   = std::numeric_limits<uint32_t>::max();
	static constexpr auto INVALID_DIST = std::numeric_limits<float>::infinity();

	MeshInstance3D *mesh_instance = nullptr;
	uint32_t surface_id           = INVALID_ID;
	uint32_t index_id             = INVALID_ID;
	float ray_origin_dist         = INVALID_DIST;

	int32_t vertex_ids[3];

	GDCLASS(MeshSurfaceIndex, RefCounted);

	public:
	MeshSurfaceIndex() = default;
	MeshSurfaceIndex(MeshInstance3D *_mesh_instance, uint32_t _surface_id, uint32_t _index_id, float _ray_origin_dist);

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
};

class TriangleRaySelect : public RefCounted
{
	GDCLASS(TriangleRaySelect, RefCounted);

	struct Ray
	{
		Vector3 Origin;
		Vector3 Normal;
	};

	struct SurfaceData
	{
		uint32_t IndexCount = 0;

		uint8_t IndexStride  = 0;
		uint8_t VertexStride = 0;

		RID IndexStorageBuffer  = RID();
		RID VertexStorageBuffer = RID();
		RID SurfaceUniformSet   = RID();
	};

	struct TriangleRaySelectShader
	{
		struct SelectedVertex
		{
			uint32_t triangle_index;
			uint32_t origin_dist;

			uint32_t vertex_ids[3];
			uint32_t pad_0;

			static constexpr auto INVALID_ID = std::numeric_limits<uint32_t>::max();
			static constexpr auto MAX_DIST   = std::numeric_limits<uint32_t>::max();
		};

		struct Params
		{
			uint32_t IndexCount;
			uint32_t IndexStride;
			uint32_t VertexStride;
			uint32_t pad0;

			float RayOrigin[3];
			uint32_t pad1;

			float RayNormal[3];
			uint32_t pad2;
		};

		TriangleRaySelectShaderRD shader;
		RID version;
		RID shader_version;
		RID pipeline;

		static constexpr uint32_t NUM_COMPUTE_SHADERS = 64;
	};

	using mesh_storage_t = RendererRD::MeshStorage;

	public:
	TriangleRaySelect();
	~TriangleRaySelect();

	static void vk_extensions_request_atomic();

	static RID get_mesh_storage_instance(MeshInstance3D *mesh_instance);
	static mesh_storage_t::MeshInstance *get_mesh_vertex_data(MeshInstance3D *mesh_instance);
	static mesh_storage_t::MeshInstance *get_mesh_vertex_data(RID mesh_storage_instance_id);

	static void _bind_methods();

	Ref<MeshSurfaceIndex> select_triangle_from_meshes(const Array &mesh_instances, const Camera3D *camera,
	                                                  const Point2i &pixel);

	Ref<MeshSurfaceIndex> select_triangle_from_mesh(MeshInstance3D *mesh_instance, const Camera3D *camera,
	                                                const Point2i &pixel);

	Ref<MeshSurfaceIndex> select_triangle_from_mesh(MeshInstance3D *mesh_instance, const Vector3 &ray_origin,
	                                                const Vector3 &ray_normal);

	SurfaceData create_mesh_instance_surface_data(const MeshInstance3D &mesh_instance, size_t surface_id,
	                                              mesh_storage_t::MeshInstance *mesh_instance_data);

	PackedVector3Array get_triangle_vertices(const Ref<MeshSurfaceIndex> &mesh_surface_index);

	private:
	TriangleRaySelectShader _triangle_ray_select_shader;

	std::map<const mesh_storage_t::MeshInstance::Surface *, SurfaceData> _surface_uniform_set;
	RID _selected_vertex_buffer;
	RID _selected_vertex_uniform_set;

	static RID generate_index_array_storage_buffer(const Mesh &mesh, int surface_id);
};

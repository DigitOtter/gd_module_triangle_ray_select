#pragma once

#include "mesh_triangle_point.h"
#include "shaders/triangle_ray_select.glsl.gen.h"
#include "triangle_transform.h"

#include <core/object/ref_counted.h>
#include <scene/3d/camera_3d.h>
#include <scene/3d/mesh_instance_3d.h>
#include <servers/rendering/renderer_rd/storage_rd/mesh_storage.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <utility>

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

			float point_on_triangle[3];
			uint32_t pad_1;

			static constexpr auto INVALID_ID = std::numeric_limits<uint32_t>::max();
			static constexpr auto MAX_DIST   = std::numeric_limits<uint32_t>::max();
		};

		struct Params
		{
			uint32_t IndexCount;
			uint32_t IndexStride;
			uint32_t VertexStride;
			uint32_t pad_0;

			float RayOrigin[3];
			uint32_t pad_1;

			float RayNormal[3];
			uint32_t pad_2;
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

	static RID get_mesh_instance_storage_instance(MeshInstance3D *mesh_instance);
	static mesh_storage_t::MeshInstance *get_mesh_instance_vertex_data(MeshInstance3D *mesh_instance);
	static mesh_storage_t::MeshInstance *get_mesh_instance_vertex_data(RID mesh_storage_instance_id);

	static RID get_mesh_storage_instance(MeshInstance3D *mesh_instance);
	static mesh_storage_t::Mesh *get_mesh_vertex_data(MeshInstance3D *mesh_instance);
	static mesh_storage_t::Mesh *get_mesh_vertex_data(RID mesh_storage_instance_id);

	static void _bind_methods();

	Ref<MeshTrianglePoint> select_triangle_from_meshes(const Array &mesh_instances, const Camera3D *camera,
	                                                   const Point2i &pixel);

	Ref<MeshTrianglePoint> select_triangle_from_mesh(MeshInstance3D *mesh_instance, const Camera3D *camera,
	                                                 const Point2i &pixel);

	Ref<MeshTrianglePoint> select_triangle_from_mesh(MeshInstance3D *mesh_instance, const Vector3 &ray_origin,
	                                                 const Vector3 &ray_normal);

	SurfaceData create_mesh_instance_surface_data(const MeshInstance3D &mesh_instance, size_t surface_id,
	                                              mesh_storage_t::MeshInstance *mesh_instance_data);
	SurfaceData create_mesh_surface_data(const MeshInstance3D &mesh_instance, size_t surface_id,
	                                     mesh_storage_t::Mesh *mesh_data);

	PackedVector3Array get_triangle_vertices(const Ref<MeshTrianglePoint> &mesh_triangle_point);

	Ref<TriangleTransform> get_triangle_transform(const Ref<MeshTrianglePoint> &mesh_triangle_point,
	                                              const Transform3D &point_tf);

	Ref<TriangleTransform> get_triangle_transform(const PackedVector3Array &triangle,
	                                              const Transform3D &point_tf) const;

	private:
	TriangleRaySelectShader _triangle_ray_select_shader;

	std::map<const mesh_storage_t::MeshInstance::Surface *, SurfaceData> _instance_surface_uniform_set;
	std::map<const mesh_storage_t::Mesh::Surface *, SurfaceData> _mesh_surface_uniform_set;

	RID _selected_vertex_buffer;
	RID _selected_vertex_uniform_set;

	static RID generate_index_array_storage_buffer(const Mesh &mesh, int surface_id);
	static std::pair<RID, uint8_t> generate_vertex_array_storage_buffer(const Mesh &mesh, int surface_id);
};

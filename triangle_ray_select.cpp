#include "triangle_ray_select.h"

#include "register_vk_extensions.h"

#include <servers/rendering/renderer_geometry_instance.h>
#include <servers/rendering/renderer_scene_cull.h>
#include <servers/rendering/rendering_device.h>
#include <servers/rendering/rendering_server_globals.h>
#include <servers/rendering_server.h>

#include <assert.h>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <strings.h>

TriangleRaySelect::TriangleRaySelect()
{
	RenderingDevice *const prd = RD::get_singleton();

	// Setup triangle_select shader
	Vector<String> shader_defines;
	shader_defines.push_back("");
	this->_triangle_ray_select_shader.shader.initialize(shader_defines);
	this->_triangle_ray_select_shader.version = this->_triangle_ray_select_shader.shader.version_create();
	this->_triangle_ray_select_shader.shader_version =
		this->_triangle_ray_select_shader.shader.version_get_shader(this->_triangle_ray_select_shader.version, 0);
	this->_triangle_ray_select_shader.pipeline =
		prd->compute_pipeline_create(this->_triangle_ray_select_shader.shader_version);

	// Create SelectedVertex buffer
	Vector<uint8_t> data;
	data.resize(sizeof(TriangleRaySelectShader::SelectedVertex));
	data.fill(0);
	this->_selected_vertex_buffer = prd->storage_buffer_create(sizeof(TriangleRaySelectShader::SelectedVertex), data);

	// Create SelectedVertex uniform set
	Vector<RD::Uniform> uniforms;
	{
		RD::Uniform u;
		u.binding      = 0;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(this->_selected_vertex_buffer);
		uniforms.push_back(u);
	}
	this->_selected_vertex_uniform_set =
		RD::get_singleton()->uniform_set_create(uniforms, this->_triangle_ray_select_shader.shader_version, 1);
}

TriangleRaySelect::~TriangleRaySelect()
{
	RenderingDevice *const prd = RD::get_singleton();
	if(this->_triangle_ray_select_shader.pipeline.is_valid())
	{
		prd->free(this->_triangle_ray_select_shader.pipeline);
		this->_triangle_ray_select_shader.pipeline = RID();
	}

	if(this->_triangle_ray_select_shader.version.is_valid())
	{
		this->_triangle_ray_select_shader.shader.version_free(this->_triangle_ray_select_shader.version);
		this->_triangle_ray_select_shader.shader_version = RID();
		this->_triangle_ray_select_shader.version        = RID();
	}

	for(auto &surface_rid: this->_instance_surface_uniform_set)
	{
		prd->free(surface_rid.second.surface_uniform_set);
		prd->free(surface_rid.second.index_storage_buffer);
	}
	this->_instance_surface_uniform_set.clear();

	for(auto &surface_rid: this->_mesh_surface_uniform_set)
	{
		prd->free(surface_rid.second.surface_uniform_set);
		prd->free(surface_rid.second.index_storage_buffer);
		prd->free(surface_rid.second.vertex_storage_buffer);
	}
	this->_mesh_surface_uniform_set.clear();

	if(this->_selected_vertex_uniform_set.is_valid())
	{
		prd->free(this->_selected_vertex_uniform_set);
		this->_selected_vertex_uniform_set = RID();
	}

	if(this->_selected_vertex_buffer.is_valid())
	{
		prd->free(this->_selected_vertex_buffer);
		this->_selected_vertex_buffer = RID();
	}
}

void TriangleRaySelect::vk_extensions_request_atomic()
{
	// Enables Vulkan's VK_KHR_SHADER_ATOMIC_INT64 extension, required for comparing distances of triangles along a
	// ray to the ray's origin.
	// Due to problems with redefinition of Window, the function is in a separate cpp file
	::vk_extensions_request_atomic();
}

void TriangleRaySelect::_bind_methods()
{
	{
		using sel_tri_fcn_t =
			Ref<MeshTrianglePoint> (TriangleRaySelect::*)(const Array &, const Camera3D *, const Point2i &);
		ClassDB::bind_method(D_METHOD("select_triangle_from_meshes_cam", "mesh_instances_array", "camera", "pixel"),
		                     (sel_tri_fcn_t)&TriangleRaySelect::select_triangle_from_meshes);
	}

	{
		using sel_tri_fcn_t =
			Ref<MeshTrianglePoint> (TriangleRaySelect::*)(const Array &, const Vector3 &, const Vector3 &);
		ClassDB::bind_method(
			D_METHOD("select_triangle_from_meshes", "mesh_instances_array", "ray_origin", "ray_normal"),
			(sel_tri_fcn_t)&TriangleRaySelect::select_triangle_from_meshes);
	}

	{
		using sel_tri_fcn_t =
			Ref<MeshTrianglePoint> (TriangleRaySelect::*)(MeshInstance3D *, const Camera3D *, const Point2i &);
		ClassDB::bind_method(D_METHOD("select_triangle_from_mesh_cam", "mesh_instance", "camera", "pixel"),
		                     (sel_tri_fcn_t)&TriangleRaySelect::select_triangle_from_mesh);
	}

	{
		using sel_tri_fcn_t =
			Ref<MeshTrianglePoint> (TriangleRaySelect::*)(MeshInstance3D *, const Vector3 &, const Vector3 &);
		ClassDB::bind_method(D_METHOD("select_triangle_from_mesh", "mesh_instance", "ray_origin", "ray_normal"),
		                     (sel_tri_fcn_t)&TriangleRaySelect::select_triangle_from_mesh);
	}

	ClassDB::bind_method(D_METHOD("get_triangle_vertices", "mesh_triangle_point"),
	                     &TriangleRaySelect::get_triangle_vertices);

	{
		using get_triangle_transform_t =
			Ref<TriangleTransform> (TriangleRaySelect::*)(const Ref<MeshTrianglePoint> &, const Transform3D &);
		ClassDB::bind_method(D_METHOD("get_triangle_transform_msi", "mesh_triangle_point", "point"),
		                     (get_triangle_transform_t)&TriangleRaySelect::get_triangle_transform);
	}

	{
		using get_triangle_transform_t =
			Ref<TriangleTransform> (TriangleRaySelect::*)(const PackedVector3Array &triangle, const Transform3D &point)
				const;
		ClassDB::bind_method(D_METHOD("get_triangle_transform", "triangle", "point"),
		                     (get_triangle_transform_t)&TriangleRaySelect::get_triangle_transform);
	}
}

Ref<MeshTrianglePoint> TriangleRaySelect::select_triangle_from_meshes(const Array &mesh_instances,
                                                                      const Camera3D *camera, const Point2i &pixel)
{
	return this->select_triangle_from_meshes(mesh_instances, camera->project_ray_origin(pixel),
	                                         camera->project_local_ray_normal(pixel));
}

Ref<MeshTrianglePoint> TriangleRaySelect::select_triangle_from_meshes(const Array &mesh_instances,
                                                                      const Vector3 &ray_origin,
                                                                      const Vector3 &ray_normal)
{
	Ref<MeshTrianglePoint> ret(memnew(MeshTrianglePoint()));
	const size_t num_meshes = mesh_instances.size();
	for(size_t i = 0; i < num_meshes; ++i)
	{
		Object *o = mesh_instances.get(i).get_validated_object();
		ERR_CONTINUE_MSG(o == nullptr || o->get_class() != MeshInstance3D::get_class_static(),
		                 "Array contains an item that is not a MeshInstance3D");

		MeshInstance3D *mesh_instance = static_cast<MeshInstance3D *>(o);
		const Transform3D &global_tf  = mesh_instance->get_global_transform();

		Ref<MeshTrianglePoint> cur_res = this->select_triangle_from_mesh(mesh_instance, global_tf.xform_inv(ray_origin),
		                                                                 global_tf.basis.xform_inv(ray_normal));

		if(cur_res->ray_origin_dist < ret->ray_origin_dist)
		{
			ret = cur_res;
		}
	}

	return ret;
}

Ref<MeshTrianglePoint> TriangleRaySelect::select_triangle_from_mesh(MeshInstance3D *mesh_instance,
                                                                    const Camera3D *camera, const Point2i &pixel)
{
	return this->select_triangle_from_mesh(mesh_instance, camera->project_ray_origin(pixel),
	                                       camera->project_ray_normal(pixel));
}

Ref<MeshTrianglePoint> TriangleRaySelect::select_triangle_from_mesh(MeshInstance3D *mesh_instance,
                                                                    const Vector3 &ray_origin,
                                                                    const Vector3 &ray_normal)
{
	ERR_FAIL_COND_V(mesh_instance == nullptr, Ref(new MeshTrianglePoint()));

	RID mesh_rid = get_mesh_storage_instance(mesh_instance);
	ERR_FAIL_COND_V(mesh_rid.is_null(), Ref(new MeshTrianglePoint()));

	// Compute shader defaults
	constexpr TriangleRaySelectShader::SelectedVertex default_selected_vertex{
		TriangleRaySelectShader::SelectedVertex::MAX_DIST};

	// Final result
	constexpr auto MAX_DIST = std::numeric_limits<uint32_t>::max();
	Ref<MeshTrianglePoint> ret(memnew(MeshTrianglePoint));
	ret->mesh_instance       = mesh_instance;
	uint32_t ray_origin_dist = MAX_DIST;

	// Ray parameters
	TriangleRaySelectShader::Params params{};
	params.ray_origin[0] = ray_origin[0];
	params.ray_origin[1] = ray_origin[1];
	params.ray_origin[2] = ray_origin[2];

	params.ray_normal[0] = ray_normal[0];
	params.ray_normal[1] = ray_normal[1];
	params.ray_normal[2] = ray_normal[2];

	RenderingDevice *const prd   = RD::get_singleton();
	mesh_storage_t &mesh_storage = *mesh_storage_t::get_singleton();

	// Mesh vertex data is stored at one of two locations. If a mesh contains blend shapes or needs a skeleton, the
	// vertex buffer is stored in a mesh_storage_t::MeshInstance struct (owned by mesh_storage_t::mesh_instance_owner).
	// If not, it's stored in a mesh_storage_t::Mesh struct (owned by mesh_storage_t::mesh_owner)
	const bool needs_instance =
		mesh_storage.mesh_needs_instance(mesh_rid, !mesh_instance->get_skeleton_path().is_empty());

	const size_t surface_count = mesh_instance->get_surface_override_material_count();
	for(size_t i = 0; i < surface_count; ++i)
	{
		// Get or create surface uniform set containing index and vertex buffer
		const SurfaceData *psurface_data;
		if(needs_instance)
		{
			mesh_storage_t::MeshInstance *pmesh_instance_data = get_mesh_instance_vertex_data(mesh_instance);
			if(auto uniform_it = this->_instance_surface_uniform_set.find(pmesh_instance_data->surfaces.ptr() + i);
			   uniform_it != this->_instance_surface_uniform_set.end())
			{
				psurface_data = &uniform_it->second;
			}
			else
			{
				SurfaceData new_surface_data =
					this->create_mesh_instance_surface_data(*mesh_instance, i, pmesh_instance_data);
				psurface_data = &this->_instance_surface_uniform_set
				                     .emplace(pmesh_instance_data->surfaces.ptr() + i, new_surface_data)
				                     .first->second;
			}
		}
		else
		{
			mesh_storage_t::Mesh *pmesh_instance_data = get_mesh_vertex_data(mesh_instance);
			if(auto uniform_it = this->_mesh_surface_uniform_set.find(pmesh_instance_data->surfaces[i]);
			   uniform_it != this->_mesh_surface_uniform_set.end())
			{
				psurface_data = &uniform_it->second;
			}
			else
			{
				SurfaceData new_surface_data = this->create_mesh_surface_data(*mesh_instance, i, pmesh_instance_data);
				psurface_data =
					&this->_mesh_surface_uniform_set.emplace(pmesh_instance_data->surfaces[i], new_surface_data)
						 .first->second;
			}
		}

		// Reset selected_vertex_buffer to default
		prd->buffer_update(this->_selected_vertex_buffer, 0, sizeof(TriangleRaySelectShader::SelectedVertex),
		                   &default_selected_vertex);

		// Prepare compute shader
		RD::ComputeListID compute_list_id = RD::get_singleton()->compute_list_begin();

		prd->compute_list_bind_uniform_set(compute_list_id, psurface_data->surface_uniform_set, 0);
		prd->compute_list_bind_uniform_set(compute_list_id, this->_selected_vertex_uniform_set, 1);

		prd->compute_list_bind_compute_pipeline(compute_list_id, this->_triangle_ray_select_shader.pipeline);

		// Set compute shader params
		params.index_count   = psurface_data->index_count;
		params.index_stride  = psurface_data->index_stride;
		params.vertex_stride = psurface_data->vertex_stride;

		prd->compute_list_set_push_constant(compute_list_id, &params, sizeof(TriangleRaySelectShader::Params));

		// Run compute shader
		const uint32_t num_triangles = params.index_stride == 3 ? params.index_count / 3 : params.index_count / 2 + 1;
		prd->compute_list_dispatch_threads(compute_list_id, num_triangles, 1, 1);

		// Wait for compute shader finish
		prd->compute_list_end();

		// Retrieve data
		const Vector<uint8_t> buf_dat =
			prd->buffer_get_data(this->_selected_vertex_buffer, 0, sizeof(TriangleRaySelectShader::SelectedVertex));

		const auto *psel_vertex = (const TriangleRaySelectShader::SelectedVertex *)buf_dat.ptr();
		if(psel_vertex->origin_dist < ray_origin_dist)
		{
			ray_origin_dist = psel_vertex->origin_dist;
			ret->surface_id = i;
			memcpy(ret->vertex_ids, psel_vertex->vertex_ids, sizeof(float) * 3);
			ret->point_on_triangle[0] = psel_vertex->point_on_triangle[0];
			ret->point_on_triangle[1] = psel_vertex->point_on_triangle[1];
			ret->point_on_triangle[2] = psel_vertex->point_on_triangle[2];
		}
	}

	if(ray_origin_dist == MAX_DIST)
	{
		ret->ray_origin_dist = MeshTrianglePoint::INVALID_DIST;
	}
	else
	{
		ret->ray_origin_dist = std::sqrt(float(ray_origin_dist) / 1000.0f);
	}
	return ret;
}

TriangleRaySelect::SurfaceData TriangleRaySelect::create_mesh_instance_surface_data(
	const MeshInstance3D &mesh_instance, size_t surface_id, mesh_storage_t::MeshInstance *mesh_instance_data) const
{
	assert(surface_id < mesh_instance_data->surfaces.size());
	const mesh_storage_t::MeshInstance::Surface *const psurface = mesh_instance_data->surfaces.ptr() + surface_id;
	const mesh_storage_t::Mesh::Surface *const pmesh_surface    = mesh_instance_data->mesh->surfaces[surface_id];

	SurfaceData surface_data;
	surface_data.index_storage_buffer =
		generate_index_array_storage_buffer(*mesh_instance.get_mesh().ptr(), surface_id);

	// From mesh_storage.cpp
	const bool has_normal               = pmesh_surface->format & RS::ARRAY_FORMAT_NORMAL;
	const bool has_tangent              = pmesh_surface->format & RS::ARRAY_FORMAT_TANGENT;
	const uint8_t normal_tangent_stride = (has_normal ? 1 : 0) + (has_tangent ? 1 : 0);

	surface_data.vertex_storage_buffer = psurface->vertex_buffer[psurface->current_buffer];
	surface_data.vertex_stride =
		(pmesh_surface->vertex_buffer_size / pmesh_surface->vertex_count) / 4 - normal_tangent_stride;

	Vector<RD::Uniform> uniforms;
	{
		RD::Uniform u;
		u.binding      = 0;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.index_storage_buffer);
		uniforms.push_back(u);
	}

	{
		RD::Uniform u;
		u.binding      = 1;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.vertex_storage_buffer);
		uniforms.push_back(u);
	}

	surface_data.index_count = pmesh_surface->index_count;

	assert(pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ||
	       pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLE_STRIP);
	surface_data.index_stride = pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ? 3 : 2;

	surface_data.surface_uniform_set =
		RD::get_singleton()->uniform_set_create(uniforms, this->_triangle_ray_select_shader.shader_version, 0);

	return surface_data;
}

TriangleRaySelect::SurfaceData TriangleRaySelect::create_mesh_surface_data(const MeshInstance3D &mesh_instance,
                                                                           size_t surface_id,
                                                                           mesh_storage_t::Mesh *mesh_data) const
{
	assert(surface_id < mesh_data->surface_count);
	const mesh_storage_t::Mesh::Surface *const pmesh_surface = mesh_data->surfaces[surface_id];

	SurfaceData surface_data;
	surface_data.index_storage_buffer =
		generate_index_array_storage_buffer(*mesh_instance.get_mesh().ptr(), surface_id);
	std::tie(surface_data.vertex_storage_buffer, surface_data.vertex_stride) =
		generate_vertex_array_storage_buffer(*mesh_instance.get_mesh().ptr(), surface_id);

	Vector<RD::Uniform> uniforms;
	{
		RD::Uniform u;
		u.binding      = 0;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.index_storage_buffer);
		uniforms.push_back(u);
	}

	{
		RD::Uniform u;
		u.binding      = 1;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.vertex_storage_buffer);
		uniforms.push_back(u);
	}

	surface_data.index_count = pmesh_surface->index_count;

	assert(pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ||
	       pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLE_STRIP);
	surface_data.index_stride = pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ? 3 : 2;

	surface_data.surface_uniform_set =
		RD::get_singleton()->uniform_set_create(uniforms, this->_triangle_ray_select_shader.shader_version, 0);

	return surface_data;
}

PackedVector3Array TriangleRaySelect::get_triangle_vertices(const Ref<MeshTrianglePoint> &mesh_triangle_point)
{
	if(unlikely(mesh_triangle_point.is_null() || mesh_triangle_point->mesh_instance == nullptr))
		return PackedVector3Array();

	const SurfaceData *psurf_data;
	if(mesh_storage_t::get_singleton()->mesh_needs_instance(
		   mesh_triangle_point->mesh_instance->get_mesh()->get_rid(),
		   !mesh_triangle_point->mesh_instance->get_skeleton_path().is_empty()))
	{
		const mesh_storage_t::MeshInstance *pmesh_data =
			TriangleRaySelect::get_mesh_instance_vertex_data(mesh_triangle_point->mesh_instance);
		ERR_FAIL_COND_V(!pmesh_data, PackedVector3Array());

		const mesh_storage_t::MeshInstance::Surface *psurface = &pmesh_data->surfaces[mesh_triangle_point->surface_id];

		const auto surf_dat_it = this->_instance_surface_uniform_set.find(psurface);
		ERR_FAIL_COND_V(surf_dat_it == this->_instance_surface_uniform_set.end(), PackedVector3Array());

		psurf_data = &surf_dat_it->second;
	}
	else
	{
		const mesh_storage_t::Mesh *pmesh_data =
			TriangleRaySelect::get_mesh_vertex_data(mesh_triangle_point->mesh_instance);
		ERR_FAIL_COND_V(!pmesh_data, PackedVector3Array());

		const mesh_storage_t::Mesh::Surface *psurface = pmesh_data->surfaces[mesh_triangle_point->surface_id];

		const auto surf_dat_it = this->_mesh_surface_uniform_set.find(psurface);
		ERR_FAIL_COND_V(surf_dat_it == this->_mesh_surface_uniform_set.end(), PackedVector3Array());

		psurf_data = &surf_dat_it->second;
	}

	RD *const prd = RD::get_singleton();

	PackedVector3Array ret;
	ret.resize(3);
	for(size_t i = 0; i < 3; ++i)
	{
		const int32_t vector_index = mesh_triangle_point->vertex_ids[i];
		const Vector<uint8_t> vector_buffer =
			prd->buffer_get_data(psurf_data->vertex_storage_buffer,
		                         vector_index * psurf_data->vertex_stride * sizeof(float), 3 * sizeof(float));
		ret.write[i][0] = *(const float *)(vector_buffer.ptr() + 0 * sizeof(float));
		ret.write[i][1] = *(const float *)(vector_buffer.ptr() + 1 * sizeof(float));
		ret.write[i][2] = *(const float *)(vector_buffer.ptr() + 2 * sizeof(float));
	}

	return ret;
}

Ref<TriangleTransform> TriangleRaySelect::get_triangle_transform(const Ref<MeshTrianglePoint> &mesh_triangle_point,
                                                                 const Transform3D &point_tf)
{
	const PackedVector3Array &triangle = this->get_triangle_vertices(mesh_triangle_point);
	return get_triangle_transform(triangle, point_tf);
}

Ref<TriangleTransform> TriangleRaySelect::get_triangle_transform(const PackedVector3Array &triangle,
                                                                 const Transform3D &point_tf) const
{
	return Ref(memnew(TriangleTransform(triangle, point_tf)));
}

RID TriangleRaySelect::generate_index_array_storage_buffer(const Mesh &mesh, int surface_id)
{
	Array surface_arrays    = mesh.surface_get_arrays(surface_id);
	Vector<int> index_array = surface_arrays[RS::ARRAY_INDEX];
	assert(index_array.size() > 0);

	// Create storage buffer
	RenderingDevice *const prd = RD::get_singleton();
	RID storage_buffer = prd->storage_buffer_create(sizeof(int) * index_array.size(), index_array.to_byte_array());
	return storage_buffer;
}

std::pair<RID, uint8_t> TriangleRaySelect::generate_vertex_array_storage_buffer(const Mesh &mesh, int surface_id)
{
	Array surface_arrays            = mesh.surface_get_arrays(surface_id);
	PackedVector3Array vertex_array = surface_arrays[RS::ARRAY_VERTEX];
	assert(vertex_array.size() > 0);

	// Create storage buffer
	RenderingDevice *const prd = RD::get_singleton();
	const auto vert_byte_array = vertex_array.to_byte_array();
	RID storage_buffer         = prd->storage_buffer_create(vert_byte_array.size(), vert_byte_array);

	const uint8_t vertex_stride = (vert_byte_array.size() / vertex_array.size()) / 4;

	return std::make_pair(storage_buffer, vertex_stride);
}

RID TriangleRaySelect::get_mesh_instance_storage_instance(MeshInstance3D *mesh_instance)
{
	// For future reference:
	// The following describes how to get an up-to-date vertex buffer from a mesh_instance pointer.
	// The MeshInstance3D RID is owned by RenderingMethod (RSG::scene). Currently, RendererSceneCull is the derived
	// class with the actual implementations.
	// MeshInstance3D's RID data is RendererSceneCull::Instance, which we can get from
	// RendererSceneCull::instance_owner. From there, we can access MeshInstance3D's base data (Instance::base_data).
	// The base_data is formatted as RendererSceneCull::InstanceGeometryData, which gets us our desired
	// InstanceGeometryData::mesh_instance.

	// Get RendererSceneCull
	auto *const pscene_cull = dynamic_cast<RendererSceneCull *>(RSG::scene);
	ERR_FAIL_NULL_V(pscene_cull, RID());
	assert(pscene_cull->instance_owner.owns(mesh_instance->get_instance()));

	// Get mesh_instance RID data
	RendererSceneCull::Instance *instance = pscene_cull->instance_owner.get_or_null(mesh_instance->get_instance());
	ERR_FAIL_NULL_V(instance, RID());

	// Get MeshInstance3D base_data from instance
	ERR_FAIL_COND_V(instance->base_type != RS::InstanceType::INSTANCE_MESH, RID());
	auto *const geom = static_cast<RendererSceneCull::InstanceGeometryData *>(instance->base_data);

	// Get mesh_id from base_data
	return dynamic_cast<RenderGeometryInstanceBase *>(geom->geometry_instance)->mesh_instance;
}

TriangleRaySelect::mesh_storage_t::MeshInstance *TriangleRaySelect::get_mesh_instance_vertex_data(
	MeshInstance3D *mesh_instance)
{
	RID mesh_storage_instance_id = get_mesh_instance_storage_instance(mesh_instance);
	ERR_FAIL_COND_V(!mesh_storage_instance_id.is_valid(), nullptr);

	return get_mesh_instance_vertex_data(mesh_storage_instance_id);
}

TriangleRaySelect::mesh_storage_t::MeshInstance *TriangleRaySelect::get_mesh_instance_vertex_data(
	RID mesh_storage_instance_id)
{
	// Check that mesh_id is really owned by mesh_storage
	mesh_storage_t &mesh_storage = *mesh_storage_t::get_singleton();
	auto &mesh_instance_owner    = mesh_storage.mesh_instance_owner;
	assert(mesh_instance_owner.owns(mesh_storage_instance_id));

	return mesh_instance_owner.get_or_null(mesh_storage_instance_id);
}

RID TriangleRaySelect::get_mesh_storage_instance(MeshInstance3D *mesh_instance)
{
	assert(mesh_storage_t::get_singleton()->mesh_owner.owns(mesh_instance->get_mesh()->get_rid()));
	return mesh_instance->get_mesh()->get_rid();
}

TriangleRaySelect::mesh_storage_t::Mesh *TriangleRaySelect::get_mesh_vertex_data(MeshInstance3D *mesh_instance)
{
	RID mesh_storage_instance_id = get_mesh_storage_instance(mesh_instance);
	ERR_FAIL_COND_V(!mesh_storage_instance_id.is_valid(), nullptr);

	return get_mesh_vertex_data(mesh_storage_instance_id);
}

TriangleRaySelect::mesh_storage_t::Mesh *TriangleRaySelect::get_mesh_vertex_data(RID mesh_storage_instance_id)
{
	// Check that mesh_id is really owned by mesh_storage
	mesh_storage_t &mesh_storage = *mesh_storage_t::get_singleton();
	auto &mesh_owner             = mesh_storage.mesh_owner;
	assert(mesh_owner.owns(mesh_storage_instance_id));

	return mesh_owner.get_or_null(mesh_storage_instance_id);
}

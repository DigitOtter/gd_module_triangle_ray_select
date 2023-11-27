#include "triangle_ray_select.h"

#include "core/error/error_macros.h"
#include "core/math/vector3.h"
#include "core/templates/vector.h"
#include "core/variant/array.h"
#include "core/variant/variant.h"
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

void MeshSurfaceIndex::_bind_methods()
{
	ClassDB::bind_method(D_METHOD("set_mesh_instance", "mesh_instance"), &MeshSurfaceIndex::set_mesh_instance);
	ClassDB::bind_method(D_METHOD("get_mesh_instance"), &MeshSurfaceIndex::get_mesh_instance);
	ClassDB::add_property("MeshSurfaceIndex",
	                      PropertyInfo(Variant::Type::OBJECT, "mesh_instance", PROPERTY_HINT_OBJECT_ID),
	                      "set_mesh_instance", "get_mesh_instance");

	ClassDB::bind_method(D_METHOD("set_surface_id", "surface_id"), &MeshSurfaceIndex::set_surface_id);
	ClassDB::bind_method(D_METHOD("get_surface_id"), &MeshSurfaceIndex::get_surface_id);
	ClassDB::add_property("MeshSurfaceIndex", PropertyInfo(Variant::Type::INT, "surface_id"), "set_surface_id",
	                      "get_surface_id");

	ClassDB::bind_method(D_METHOD("set_index_id", "index_id"), &MeshSurfaceIndex::set_index_id);
	ClassDB::bind_method(D_METHOD("get_index_id"), &MeshSurfaceIndex::get_index_id);
	ClassDB::add_property("MeshSurfaceIndex", PropertyInfo(Variant::Type::INT, "index_id"), "set_index_id",
	                      "get_index_id");

	ClassDB::bind_method(D_METHOD("set_ray_origin_dist", "ray_origin_dist"), &MeshSurfaceIndex::set_ray_origin_dist);
	ClassDB::bind_method(D_METHOD("get_ray_origin_dist"), &MeshSurfaceIndex::get_ray_origin_dist);
	ClassDB::add_property("MeshSurfaceIndex", PropertyInfo(Variant::Type::FLOAT, "ray_origin_dist"),
	                      "set_ray_origin_dist", "get_ray_origin_dist");

	ClassDB::bind_method(D_METHOD("set_vertex_ids", "vertex_ids"), &MeshSurfaceIndex::set_vertex_ids);
	ClassDB::bind_method(D_METHOD("get_vertex_ids"), &MeshSurfaceIndex::get_vertex_ids);
	ClassDB::add_property("MeshSurfaceIndex", PropertyInfo(Variant::Type::PACKED_INT32_ARRAY, "vertex_ids"),
	                      "set_vertex_ids", "get_vertex_ids");
}

MeshSurfaceIndex::MeshSurfaceIndex(MeshInstance3D *_mesh_instance, uint32_t _surface_id, uint32_t _index_id,
                                   float _ray_origin_dist)
	: mesh_instance(_mesh_instance),
	  surface_id(_surface_id),
	  index_id(_index_id),
	  ray_origin_dist(_ray_origin_dist)
{}

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
		prd->free(surface_rid.second.SurfaceUniformSet);
		prd->free(surface_rid.second.IndexStorageBuffer);
	}
	this->_instance_surface_uniform_set.clear();

	for(auto &surface_rid: this->_mesh_surface_uniform_set)
	{
		prd->free(surface_rid.second.SurfaceUniformSet);
		prd->free(surface_rid.second.IndexStorageBuffer);
		prd->free(surface_rid.second.VertexStorageBuffer);
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
	::vk_extensions_request_atomic();
}

RID TriangleRaySelect::get_mesh_instance_storage_instance(MeshInstance3D *mesh_instance)
{
	// For future reference:
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

void TriangleRaySelect::_bind_methods()
{
	using sel_tri_fcn_t =
		Ref<MeshSurfaceIndex> (TriangleRaySelect::*)(MeshInstance3D *, const Camera3D *, const Point2i &);
	ClassDB::bind_method(D_METHOD("select_triangle_from_mesh", "mesh_instance", "camera", "pixel"),
	                     (sel_tri_fcn_t)&TriangleRaySelect::select_triangle_from_mesh);

	ClassDB::bind_method(D_METHOD("select_triangle_from_meshes", "mesh_instances_array", "camera", "pixel"),
	                     &TriangleRaySelect::select_triangle_from_meshes);

	ClassDB::bind_method(D_METHOD("get_triangle_vertices", "mesh_surface_index"),
	                     &TriangleRaySelect::get_triangle_vertices);
}

Ref<MeshSurfaceIndex> TriangleRaySelect::select_triangle_from_meshes(const Array &mesh_instances,
                                                                     const Camera3D *camera, const Point2i &pixel)
{
	const Vector3 ray_origin = camera->project_ray_origin(pixel);
	const Vector3 ray_normal = camera->project_ray_normal(pixel);

	Ref<MeshSurfaceIndex> ret(memnew(MeshSurfaceIndex()));
	const size_t num_meshes = mesh_instances.size();
	for(size_t i = 0; i < num_meshes; ++i)
	{
		Object *o = mesh_instances.get(i).get_validated_object();
		ERR_CONTINUE_MSG(o == nullptr || o->get_class() != MeshInstance3D::get_class_static(),
		                 "Array contains an item that is not a MeshInstance3D");

		MeshInstance3D *mesh_instance = static_cast<MeshInstance3D *>(o);
		const Transform3D &global_tf  = mesh_instance->get_global_transform();

		Ref<MeshSurfaceIndex> cur_res = this->select_triangle_from_mesh(mesh_instance, global_tf.xform_inv(ray_origin),
		                                                                global_tf.basis.xform_inv(ray_normal));

		if(cur_res->ray_origin_dist < ret->ray_origin_dist)
		{
			ret = cur_res;
		}
	}

	return ret;
}

Ref<MeshSurfaceIndex> TriangleRaySelect::select_triangle_from_mesh(MeshInstance3D *mesh_instance,
                                                                   const Camera3D *camera, const Point2i &pixel)
{
	return this->select_triangle_from_mesh(mesh_instance, camera->project_ray_origin(pixel),
	                                       camera->project_ray_normal(pixel));
}

Ref<MeshSurfaceIndex> TriangleRaySelect::select_triangle_from_mesh(MeshInstance3D *mesh_instance,
                                                                   const Vector3 &ray_origin, const Vector3 &ray_normal)
{
	ERR_FAIL_COND_V(mesh_instance == nullptr, Ref(new MeshSurfaceIndex()));

	RID mesh_rid = get_mesh_storage_instance(mesh_instance);
	ERR_FAIL_COND_V(mesh_rid.is_null(), Ref(new MeshSurfaceIndex()));

	// RID mesh_storage_instance_id                            = get_mesh_instance_storage_instance(mesh_instance);
	// mesh_storage_t::MeshInstance *const pmesh_instance_data =
	// get_mesh_instance_vertex_data(mesh_storage_instance_id);

	// Apply skeleton and blendshape deforms
	// mesh_storage_t *pmesh_storage = mesh_storage_t::get_singleton();
	// pmesh_storage->mesh_instance_check_for_update(mesh_storage_instance_id);
	// pmesh_storage->update_mesh_instances();

	// Compute shader defaults
	constexpr TriangleRaySelectShader::SelectedVertex default_selected_vertex{
		TriangleRaySelectShader::SelectedVertex::INVALID_ID, TriangleRaySelectShader::SelectedVertex::MAX_DIST};

	// Final result
	constexpr auto MAX_DIST = std::numeric_limits<uint32_t>::max();
	Ref<MeshSurfaceIndex> ret(memnew(MeshSurfaceIndex));
	ret->mesh_instance       = mesh_instance;
	uint32_t ray_origin_dist = MAX_DIST;

	// Ray parameters
	TriangleRaySelectShader::Params params{};
	params.RayOrigin[0] = ray_origin[0];
	params.RayOrigin[1] = ray_origin[1];
	params.RayOrigin[2] = ray_origin[2];

	params.RayNormal[0] = ray_normal[0];
	params.RayNormal[1] = ray_normal[1];
	params.RayNormal[2] = ray_normal[2];

	RenderingDevice *const prd   = RD::get_singleton();
	mesh_storage_t &mesh_storage = *mesh_storage_t::get_singleton();

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

		prd->compute_list_bind_uniform_set(compute_list_id, psurface_data->SurfaceUniformSet, 0);
		prd->compute_list_bind_uniform_set(compute_list_id, this->_selected_vertex_uniform_set, 1);

		prd->compute_list_bind_compute_pipeline(compute_list_id, this->_triangle_ray_select_shader.pipeline);

		// Set compute shader params
		params.IndexCount   = psurface_data->IndexCount;
		params.IndexStride  = psurface_data->IndexStride;
		params.VertexStride = psurface_data->VertexStride;

		prd->compute_list_set_push_constant(compute_list_id, &params, sizeof(TriangleRaySelectShader::Params));

		// Run compute shader
		const uint32_t num_triangles = params.IndexStride == 3 ? params.IndexCount / 3 : params.IndexCount / 2 + 1;
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
			ret->index_id   = psel_vertex->triangle_index;
			ret->surface_id = i;

			memcpy(ret->vertex_ids, psel_vertex->vertex_ids, sizeof(float) * 3);
		}
	}

	if(ray_origin_dist == MAX_DIST)
	{
		ret->ray_origin_dist = MeshSurfaceIndex::INVALID_DIST;
	}
	else
	{
		ret->ray_origin_dist = std::sqrt(float(ray_origin_dist) / 1000.0f);
	}
	return ret;
}

TriangleRaySelect::SurfaceData TriangleRaySelect::create_mesh_instance_surface_data(
	const MeshInstance3D &mesh_instance, size_t surface_id, mesh_storage_t::MeshInstance *mesh_instance_data)
{
	assert(surface_id < mesh_instance_data->surfaces.size());
	const mesh_storage_t::MeshInstance::Surface *const psurface = mesh_instance_data->surfaces.ptr() + surface_id;
	const mesh_storage_t::Mesh::Surface *const pmesh_surface    = mesh_instance_data->mesh->surfaces[surface_id];

	SurfaceData surface_data;
	surface_data.IndexStorageBuffer  = generate_index_array_storage_buffer(*mesh_instance.get_mesh().ptr(), surface_id);

	surface_data.VertexStorageBuffer = psurface->vertex_buffer;
	surface_data.VertexStride        = (pmesh_surface->vertex_buffer_size / pmesh_surface->vertex_count) / 4;

	Vector<RD::Uniform> uniforms;
	{
		RD::Uniform u;
		u.binding      = 0;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.IndexStorageBuffer);
		uniforms.push_back(u);
	}

	{
		RD::Uniform u;
		u.binding      = 1;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.VertexStorageBuffer);
		uniforms.push_back(u);
	}

	surface_data.IndexCount = pmesh_surface->index_count;

	assert(pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ||
	       pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLE_STRIP);
	surface_data.IndexStride = pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ? 3 : 2;

	surface_data.SurfaceUniformSet =
		RD::get_singleton()->uniform_set_create(uniforms, this->_triangle_ray_select_shader.shader_version, 0);

	return surface_data;
}

TriangleRaySelect::SurfaceData TriangleRaySelect::create_mesh_surface_data(const MeshInstance3D &mesh_instance,
                                                                           size_t surface_id,
                                                                           mesh_storage_t::Mesh *mesh_data)
{
	assert(surface_id < mesh_data->surface_count);
	const mesh_storage_t::Mesh::Surface *const pmesh_surface = mesh_data->surfaces[surface_id];

	SurfaceData surface_data;
	surface_data.IndexStorageBuffer = generate_index_array_storage_buffer(*mesh_instance.get_mesh().ptr(), surface_id);
	std::tie(surface_data.VertexStorageBuffer, surface_data.VertexStride) =
		generate_vertex_array_storage_buffer(*mesh_instance.get_mesh().ptr(), surface_id);

	Vector<RD::Uniform> uniforms;
	{
		RD::Uniform u;
		u.binding      = 0;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.IndexStorageBuffer);
		uniforms.push_back(u);
	}

	{
		RD::Uniform u;
		u.binding      = 1;
		u.uniform_type = RD::UNIFORM_TYPE_STORAGE_BUFFER;
		u.append_id(surface_data.VertexStorageBuffer);
		uniforms.push_back(u);
	}

	surface_data.IndexCount = pmesh_surface->index_count;

	assert(pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ||
	       pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLE_STRIP);
	surface_data.IndexStride = pmesh_surface->primitive == RenderingServer::PRIMITIVE_TRIANGLES ? 3 : 2;

	surface_data.SurfaceUniformSet =
		RD::get_singleton()->uniform_set_create(uniforms, this->_triangle_ray_select_shader.shader_version, 0);

	return surface_data;
}

PackedVector3Array TriangleRaySelect::get_triangle_vertices(const Ref<MeshSurfaceIndex> &mesh_surface_index)
{
	ERR_FAIL_COND_V(mesh_surface_index.is_null(), PackedVector3Array());
	ERR_FAIL_COND_V(mesh_surface_index->index_id == MeshSurfaceIndex::INVALID_ID, PackedVector3Array());

	const SurfaceData *psurf_data;
	if(mesh_storage_t::get_singleton()->mesh_needs_instance(
		   mesh_surface_index->mesh_instance->get_mesh()->get_rid(),
		   !mesh_surface_index->mesh_instance->get_skeleton_path().is_empty()))
	{
		const mesh_storage_t::MeshInstance *pmesh_data =
			TriangleRaySelect::get_mesh_instance_vertex_data(mesh_surface_index->mesh_instance);
		ERR_FAIL_COND_V(!pmesh_data, PackedVector3Array());

		const mesh_storage_t::MeshInstance::Surface *psurface = &pmesh_data->surfaces[mesh_surface_index->surface_id];

		const auto surf_dat_it = this->_instance_surface_uniform_set.find(psurface);
		ERR_FAIL_COND_V(surf_dat_it == this->_instance_surface_uniform_set.end(), PackedVector3Array());

		psurf_data = &surf_dat_it->second;
	}
	else
	{
		const mesh_storage_t::Mesh *pmesh_data =
			TriangleRaySelect::get_mesh_vertex_data(mesh_surface_index->mesh_instance);
		ERR_FAIL_COND_V(!pmesh_data, PackedVector3Array());

		const mesh_storage_t::Mesh::Surface *psurface = pmesh_data->surfaces[mesh_surface_index->surface_id];

		const auto surf_dat_it = this->_mesh_surface_uniform_set.find(psurface);
		ERR_FAIL_COND_V(surf_dat_it == this->_mesh_surface_uniform_set.end(), PackedVector3Array());

		psurf_data = &surf_dat_it->second;
	}

	RD *const prd = RD::get_singleton();

	PackedVector3Array ret;
	ret.resize(3);
	for(size_t i = 0; i < 3; ++i)
	{
		const int32_t vector_index = mesh_surface_index->vertex_ids[i];
		const Vector<uint8_t> vector_buffer =
			prd->buffer_get_data(psurf_data->VertexStorageBuffer,
		                         vector_index * psurf_data->VertexStride * sizeof(float), 3 * sizeof(float));
		ret.write[i][0] = *(const float *)(vector_buffer.ptr() + 0 * sizeof(float));
		ret.write[i][1] = *(const float *)(vector_buffer.ptr() + 1 * sizeof(float));
		ret.write[i][2] = *(const float *)(vector_buffer.ptr() + 2 * sizeof(float));
	}

	return ret;
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
	return std::make_pair(storage_buffer, (vert_byte_array.size() / vertex_array.size()) / 4);
}

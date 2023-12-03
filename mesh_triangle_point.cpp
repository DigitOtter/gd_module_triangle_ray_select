#include "mesh_triangle_point.h"

void MeshTrianglePoint::_bind_methods()
{
	ClassDB::bind_method(D_METHOD("set_mesh_instance", "mesh_instance"), &MeshTrianglePoint::set_mesh_instance);
	ClassDB::bind_method(D_METHOD("get_mesh_instance"), &MeshTrianglePoint::get_mesh_instance);
	ClassDB::add_property("MeshTrianglePoint",
	                      PropertyInfo(Variant::Type::OBJECT, "mesh_instance", PROPERTY_HINT_OBJECT_ID),
	                      "set_mesh_instance", "get_mesh_instance");

	ClassDB::bind_method(D_METHOD("set_surface_id", "surface_id"), &MeshTrianglePoint::set_surface_id);
	ClassDB::bind_method(D_METHOD("get_surface_id"), &MeshTrianglePoint::get_surface_id);
	ClassDB::add_property("MeshTrianglePoint", PropertyInfo(Variant::Type::INT, "surface_id"), "set_surface_id",
	                      "get_surface_id");

	ClassDB::bind_method(D_METHOD("set_index_id", "index_id"), &MeshTrianglePoint::set_index_id);
	ClassDB::bind_method(D_METHOD("get_index_id"), &MeshTrianglePoint::get_index_id);
	ClassDB::add_property("MeshTrianglePoint", PropertyInfo(Variant::Type::INT, "index_id"), "set_index_id",
	                      "get_index_id");

	ClassDB::bind_method(D_METHOD("set_ray_origin_dist", "ray_origin_dist"), &MeshTrianglePoint::set_ray_origin_dist);
	ClassDB::bind_method(D_METHOD("get_ray_origin_dist"), &MeshTrianglePoint::get_ray_origin_dist);
	ClassDB::add_property("MeshTrianglePoint", PropertyInfo(Variant::Type::FLOAT, "ray_origin_dist"),
	                      "set_ray_origin_dist", "get_ray_origin_dist");

	ClassDB::bind_method(D_METHOD("set_vertex_ids", "vertex_ids"), &MeshTrianglePoint::set_vertex_ids);
	ClassDB::bind_method(D_METHOD("get_vertex_ids"), &MeshTrianglePoint::get_vertex_ids);
	ClassDB::add_property("MeshTrianglePoint", PropertyInfo(Variant::Type::PACKED_INT32_ARRAY, "vertex_ids"),
	                      "set_vertex_ids", "get_vertex_ids");

	ClassDB::bind_method(D_METHOD("set_point_on_triangle", "point_on_triangle"),
	                     &MeshTrianglePoint::set_point_on_triangle);
	ClassDB::bind_method(D_METHOD("get_point_on_triangle"), &MeshTrianglePoint::get_point_on_triangle);
	ClassDB::add_property("MeshTrianglePoint", PropertyInfo(Variant::Type::VECTOR3, "point_on_triangle"),
	                      "set_point_on_triangle", "get_point_on_triangle");
}

MeshTrianglePoint::MeshTrianglePoint(MeshInstance3D *_mesh_instance, uint32_t _surface_id, uint32_t _index_id,
                                     float _ray_origin_dist)
	: mesh_instance(_mesh_instance),
	  surface_id(_surface_id),
	  index_id(_index_id),
	  ray_origin_dist(_ray_origin_dist)
{}
